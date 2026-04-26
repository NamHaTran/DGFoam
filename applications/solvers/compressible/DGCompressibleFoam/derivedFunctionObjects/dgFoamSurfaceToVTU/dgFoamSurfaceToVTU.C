/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2025 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgFoamSurfaceToVTU.H"

#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "IOdictionary.H"
#include "IOobject.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "dgField.H"
#include "dgThermoConservative.H"
#include "dgAffineMapping.H"
#include "dgVtkLagrangeTools.H"
#include "basisFunctions.H"
#include "boundBox.H"
#include "refCoordTransforms.H"
#include "volFields.H"

#include <cmath>
#include <unordered_map>
#include <vector>

namespace Foam
{
namespace dgFunctionObjects
{

namespace
{

constexpr unsigned char vtkTriType = 5u;
constexpr unsigned char vtkQuadType = 9u;

struct PointMergeKey
{
    long long x;
    long long y;
    long long z;

    bool operator==(const PointMergeKey& other) const
    {
        return x == other.x && y == other.y && z == other.z;
    }
};


struct PointMergeKeyHash
{
    std::size_t operator()(const PointMergeKey& key) const
    {
        std::size_t seed = std::hash<long long>{}(key.x);
        seed ^= std::hash<long long>{}(key.y) + 0x9e3779b97f4a7c15ULL
             + (seed << 6) + (seed >> 2);
        seed ^= std::hash<long long>{}(key.z) + 0x9e3779b97f4a7c15ULL
             + (seed << 6) + (seed >> 2);
        return seed;
    }
};


inline PointMergeKey makePointMergeKey(const point& pt, const scalar tol)
{
    return
    {
        std::llround(pt.x()/tol),
        std::llround(pt.y()/tol),
        std::llround(pt.z()/tol)
    };
}


inline label trianglePointCount(const label n)
{
    return (n + 1)*(n + 2)/2;
}


inline label triangleRowOffset(const label n, const label row)
{
    return row*(n + 1) - row*(row - 1)/2;
}


inline label trianglePointIndex
(
    const label n,
    const label row,
    const label col
)
{
    return triangleRowOffset(n, row) + col;
}


inline dgCellType cellTypeFromVertexCount(const label nVertices)
{
    switch (nVertices)
    {
        case 4:
            return dgCellType::TET;

        case 5:
            return dgCellType::PYRAMID;

        case 6:
            return dgCellType::PRISM;

        case 8:
            return dgCellType::HEX;

        default:
            FatalErrorInFunction
                << "Unsupported DG cell with " << nVertices
                << " vertices." << exit(FatalError);
    }

    return dgCellType::INVALID;
}


inline void computeBasisAt
(
    const dgCellType type,
    const vector& eta,
    const label pOrder,
    List<scalar>& basis
)
{
    List<scalar> dEta1;
    List<scalar> dEta2;
    List<scalar> dEta3;

    switch (type)
    {
        case dgCellType::HEX:
            computeHexBasisAndDerivatives(eta, pOrder, basis, dEta1, dEta2, dEta3);
            break;

        case dgCellType::PRISM:
            computePrismBasisAndDerivatives(eta, pOrder, basis, dEta1, dEta2, dEta3);
            break;

        case dgCellType::PYRAMID:
            computePyramidBasisAndDerivatives(eta, pOrder, basis, dEta1, dEta2, dEta3);
            break;

        case dgCellType::TET:
            computeTetBasisAndDerivatives(eta, pOrder, basis, dEta1, dEta2, dEta3);
            break;

        default:
            FatalErrorInFunction
                << "Unsupported DG cell type " << type
                << exit(FatalError);
    }
}

} // End anonymous namespace


defineTypeNameAndDebug(dgFoamSurfaceToVTU, 0);
addToRunTimeSelectionTable(functionObject, dgFoamSurfaceToVTU, dictionary);


dgFoamSurfaceToVTU::dgFoamSurfaceToVTU
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    runTime_(runTime),
    regionName_(polyMesh::defaultRegion),
    rhoName_("rho"),
    rhoUName_("rhoU"),
    EName_("E"),
    thermoName_("dgThermoConservative"),
    n_(0),
    pOrder_(-1),
    mergeSharedPoints_(false),
    pointMergeTol_(-1),
    patchNames_(),
    patchIDs_(),
    initialized_(false)
{
    read(dict);

    Info<< type() << ':' << nl
        << "    " << name << " -> surface-vtu" << nl
        << "        initialized" << endl;
}


bool dgFoamSurfaceToVTU::read(const dictionary& dict)
{
    functionObject::read(dict);

    regionName_ = dict.getOrDefault<word>("region", polyMesh::defaultRegion);
    rhoName_ = dict.getOrDefault<word>("rho", "rho");
    rhoUName_ = dict.getOrDefault<word>("rhoU", "rhoU");
    EName_ = dict.getOrDefault<word>("E", "E");
    thermoName_ = dict.getOrDefault<word>("thermo", "dgThermoConservative");

    mergeSharedPoints_ =
        dict.getOrDefault<Switch>("mergeSharedPoints", false);
    pointMergeTol_ = dict.getOrDefault<scalar>("pointMergeTol", -1);

    if (pointMergeTol_ == 0)
    {
        FatalIOErrorInFunction(dict)
            << "Expected 'pointMergeTol' to be non-zero for "
            << "dgFoamSurfaceToVTU." << exit(FatalIOError);
    }

    if (!dict.found("patches"))
    {
        FatalIOErrorInFunction(dict)
            << "Missing required entry 'patches' for "
            << "dgFoamSurfaceToVTU." << exit(FatalIOError);
    }

    patchNames_ = dict.get<wordList>("patches");

    if (!patchNames_.size())
    {
        FatalIOErrorInFunction(dict)
            << "Expected non-empty 'patches' list for dgFoamSurfaceToVTU."
            << exit(FatalIOError);
    }

    const bool hasOutputPoints = dict.found("outputPoints");
    const bool hasLegacyN = dict.found("n");

    if (hasOutputPoints)
    {
        const label outputPoints = readLabel(dict.lookup("outputPoints"));

        if (outputPoints < 2)
        {
            FatalIOErrorInFunction(dict)
                << "Expected 'outputPoints >= 2' for dgFoamSurfaceToVTU, got "
                << outputPoints << exit(FatalIOError);
        }

        n_ = outputPoints - 1;

        if (hasLegacyN)
        {
            WarningInFunction
                << "Both 'outputPoints' and legacy alias 'n' were provided for "
                << type() << ". Using 'outputPoints=" << outputPoints
                << "' and ignoring 'n'." << endl;
        }
    }
    else if (hasLegacyN)
    {
        n_ = readLabel(dict.lookup("n"));

        if (n_ < 1)
        {
            FatalIOErrorInFunction(dict)
                << "Expected legacy 'n >= 1' for dgFoamSurfaceToVTU, got "
                << n_ << exit(FatalIOError);
        }
    }
    else
    {
        n_ = 0;
    }

    patchIDs_.setSize(0);
    pOrder_ = -1;
    initialized_ = false;

    return true;
}


bool dgFoamSurfaceToVTU::execute()
{
    return true;
}


bool dgFoamSurfaceToVTU::write()
{
    initialize();

    const fileName dir = outputDir();
    mkDir(dir);

    writePiece(pieceFileName());

    if (Pstream::master() && Pstream::parRun())
    {
        writePvtu(pvtuFileName());
    }

    if (log)
    {
        Info<< type() << " wrote surface VTU samples" << endl;
    }

    return true;
}


void dgFoamSurfaceToVTU::initialize()
{
    if (initialized_)
    {
        return;
    }

    const fvMesh& meshRef = mesh();
    const dgThermoConservative& dgThermo = thermo();

    if (!dgThermo.heIsInternalEnergy())
    {
        FatalErrorInFunction
            << type()
            << " currently supports only internal-energy based thermo."
            << nl
            << "Detected energy variable: " << dgThermo.heName()
            << exit(FatalError);
    }

    IOdictionary dgSchemesDict
    (
        IOobject
        (
            "dgSchemes",
            runTime_.system(),
            meshRef,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    pOrder_ =
        readLabel(dgSchemesDict.subDict("dgDiscretization").lookup("pOrder"));

    if (n_ == 0)
    {
        const label outputPoints = max(2, pOrder_ + 1);
        n_ = outputPoints - 1;

        if (log)
        {
            Info<< type()
                << " defaulted outputPoints to " << outputPoints
                << " from pOrder=" << pOrder_ << endl;
        }
    }

    const fvBoundaryMesh& patches = meshRef.boundary();
    patchIDs_.setSize(patchNames_.size());

    forAll(patchNames_, patchNameI)
    {
        const word& targetName = patchNames_[patchNameI];
        label patchID = -1;

        forAll(patches, patchI)
        {
            if (patches[patchI].name() == targetName)
            {
                patchID = patchI;
                break;
            }
        }

        if (patchID < 0)
        {
            FatalErrorInFunction
                << "Cannot find boundary patch '" << targetName
                << "' requested by " << type() << exit(FatalError);
        }

        patchIDs_[patchNameI] = patchID;
    }

    initialized_ = true;
}


const fvMesh& dgFoamSurfaceToVTU::mesh() const
{
    if (!runTime_.foundObject<fvMesh>(regionName_))
    {
        FatalErrorInFunction
            << "Cannot find fvMesh region '" << regionName_ << "' required by "
            << type() << exit(FatalError);
    }

    return runTime_.lookupObject<fvMesh>(regionName_);
}


const dgField<scalar>& dgFoamSurfaceToVTU::rhoField() const
{
    return mesh().lookupObject<dgField<scalar>>(rhoName_);
}


const dgField<vector>& dgFoamSurfaceToVTU::rhoUField() const
{
    return mesh().lookupObject<dgField<vector>>(rhoUName_);
}


const dgField<scalar>& dgFoamSurfaceToVTU::EField() const
{
    return mesh().lookupObject<dgField<scalar>>(EName_);
}


const dgThermoConservative& dgFoamSurfaceToVTU::thermo() const
{
    return mesh().lookupObject<dgThermoConservative>(thermoName_);
}


fileName dgFoamSurfaceToVTU::outputDir() const
{
    return
        runTime_.globalPath()
      / functionObject::outputPrefix
      / name()
      / runTime_.timeName();
}


fileName dgFoamSurfaceToVTU::pieceFileName() const
{
    const fileName base = outputDir();

    if (Pstream::parRun())
    {
        return base/(name() + "_" + Foam::name(Pstream::myProcNo()) + ".vtu");
    }

    return base/(name() + ".vtu");
}


fileName dgFoamSurfaceToVTU::pvtuFileName() const
{
    return outputDir()/(name() + ".pvtu");
}


template<class Type>
Type dgFoamSurfaceToVTU::evaluateField
(
    const dgField<Type>& field,
    const label cellI,
    const List<scalar>& basis,
    const scalar highOrderScale
)
{
    const List<Type>& cellDof = field.dof()[cellI].dof();

    if (cellDof.size() != basis.size())
    {
        FatalErrorInFunction
            << "Basis size (" << basis.size()
            << ") does not match number of DoFs (" << cellDof.size()
            << ") in cell " << cellI << exit(FatalError);
    }

    Type value = pTraits<Type>::zero;

    forAll(cellDof, dofI)
    {
        const scalar modeScale = (dofI == 0 ? scalar(1) : highOrderScale);
        value += modeScale*cellDof[dofI]*basis[dofI];
    }

    return value;
}


void dgFoamSurfaceToVTU::writePiece(const fileName& filePath) const
{
    const fvMesh& meshRef = mesh();
    const dgField<scalar>& rho = rhoField();
    const dgField<vector>& rhoU = rhoUField();
    const dgField<scalar>& E = EField();
    const dgThermoConservative& thermoRef = thermo();

    const volScalarField* theta1FieldPtr =
        meshRef.foundObject<volScalarField>("theta1")
      ? &meshRef.lookupObject<volScalarField>("theta1")
      : nullptr;

    const volScalarField* theta2FieldPtr =
        meshRef.foundObject<volScalarField>("theta2")
      ? &meshRef.lookupObject<volScalarField>("theta2")
      : nullptr;

    dgVtkLagrange::GridBuffers grid;
    std::vector<scalar> rhoValues;
    std::vector<vector> rhoUValues;
    std::vector<scalar> EValues;
    std::vector<vector> UValues;
    std::vector<scalar> eValues;
    std::vector<scalar> pValues;
    std::vector<scalar> TValues;
    std::vector<scalar> aValues;
    std::vector<scalar> theta1Values;
    std::vector<scalar> theta2Values;

    scalar mergeTol = -1;
    std::unordered_map<PointMergeKey, label, PointMergeKeyHash> mergedPointMap;

    if (mergeSharedPoints_)
    {
        const scalar defaultTol =
            max(SMALL, 1e-12*boundBox(meshRef.points(), false).mag());

        mergeTol = (pointMergeTol_ > 0) ? pointMergeTol_ : defaultTol;
        mergedPointMap.reserve(1024);

        if (log)
        {
            Info<< type()
                << " merging coincident surface points with tolerance "
                << mergeTol << endl;
        }
    }

    const auto appendPoint =
    [&]
    (
        const point& physicalPoint,
        const label ownerCellI,
        const pointField& ownerCellVertices,
        const dgCellType ownerCellType,
        const scalar theta1,
        const scalar theta2
    ) -> label
    {
        if (mergeSharedPoints_)
        {
            const PointMergeKey key = makePointMergeKey(physicalPoint, mergeTol);
            const auto iter = mergedPointMap.find(key);

            if (iter != mergedPointMap.end())
            {
                return iter->second;
            }

            const label pointI = grid.points.size();
            mergedPointMap.emplace(key, pointI);
        }

        const vector ownerEta =
            mapXToEta(physicalPoint, ownerCellVertices, ownerCellType);

        List<scalar> basis;
        computeBasisAt(ownerCellType, ownerEta, pOrder_, basis);

        const scalar rhoTheta = theta1*theta2;
        const scalar rhoValue = evaluateField(rho, ownerCellI, basis, rhoTheta);
        const vector rhoUValue = evaluateField(rhoU, ownerCellI, basis, theta2);
        const scalar EValue = evaluateField(E, ownerCellI, basis, theta2);

        const scalar rhoSafe = max(rhoValue, SMALL);
        const vector UValue = rhoUValue/rhoSafe;
        const scalar eValue = EValue/rhoSafe - 0.5*magSqr(UValue);

        const label pointI = grid.points.size();
        grid.points.push_back(physicalPoint);
        rhoValues.push_back(rhoValue);
        rhoUValues.push_back(rhoUValue);
        EValues.push_back(EValue);
        UValues.push_back(UValue);
        eValues.push_back(eValue);
        TValues.push_back
        (
            thermoRef.calcTemperatureFromRhoHe(ownerCellI, rhoSafe, eValue)
        );
        pValues.push_back
        (
            thermoRef.calcPressureFromRhoHe(ownerCellI, rhoSafe, eValue)
        );
        aValues.push_back
        (
            thermoRef.calcSpeedOfSoundFromRhoHe(ownerCellI, rhoSafe, eValue)
        );
        theta1Values.push_back(theta1);
        theta2Values.push_back(theta2);

        return pointI;
    };

    const fvBoundaryMesh& patches = meshRef.boundary();

    forAll(patchIDs_, patchListI)
    {
        const label patchID = patchIDs_[patchListI];
        const fvPatch& patch = patches[patchID];
        const label start = patch.start();
        const label end = start + patch.size();

        for (label faceI = start; faceI < end; ++faceI)
        {
            const face& meshFace = meshRef.faces()[faceI];
            const label nFacePoints = meshFace.size();

            if (nFacePoints != 3 && nFacePoints != 4)
            {
                WarningInFunction
                    << "Skipping patch face " << faceI << " on patch '"
                    << patch.name() << "' because only TRI/QUAD faces are "
                    << "currently supported, got " << nFacePoints
                    << " vertices." << endl;
                continue;
            }

            pointField faceVertices(nFacePoints);

            forAll(meshFace, pointI)
            {
                faceVertices[pointI] = meshRef.points()[meshFace[pointI]];
            }

            const label ownerCellI = meshRef.faceOwner()[faceI];
            const cellShape& ownerShape = meshRef.cellShapes()[ownerCellI];
            const pointField ownerCellVertices =
                ownerShape.points(meshRef.points());
            const dgCellType ownerCellType =
                cellTypeFromVertexCount(ownerCellVertices.size());
            const scalar theta1 =
                theta1FieldPtr
              ? theta1FieldPtr->primitiveField()[ownerCellI]
              : scalar(1);
            const scalar theta2 =
                theta2FieldPtr
              ? theta2FieldPtr->primitiveField()[ownerCellI]
              : scalar(1);

            if (nFacePoints == 3)
            {
                labelList localToGlobal(trianglePointCount(n_), -1);

                for (label row = 0; row <= n_; ++row)
                {
                    for (label col = 0; col <= n_ - row; ++col)
                    {
                        const label localPointI =
                            trianglePointIndex(n_, row, col);
                        const vector2D rs
                        (
                            scalar(col)/scalar(n_),
                            scalar(row)/scalar(n_)
                        );
                        const point physicalPoint =
                            mapToPhysicalTriangle(rs, faceVertices);
                        localToGlobal[localPointI] =
                            appendPoint
                            (
                                physicalPoint,
                                ownerCellI,
                                ownerCellVertices,
                                ownerCellType,
                                theta1,
                                theta2
                            );
                    }
                }

                for (label row = 0; row < n_; ++row)
                {
                    for (label col = 0; col < n_ - row; ++col)
                    {
                        grid.connectivity.push_back
                        (
                            localToGlobal[trianglePointIndex(n_, row, col)]
                        );
                        grid.connectivity.push_back
                        (
                            localToGlobal[trianglePointIndex(n_, row, col + 1)]
                        );
                        grid.connectivity.push_back
                        (
                            localToGlobal[trianglePointIndex(n_, row + 1, col)]
                        );
                        grid.offsets.push_back(grid.connectivity.size());
                        grid.cellTypes.push_back(vtkTriType);
                        grid.partitions.push_back(Pstream::myProcNo());

                        if (row + col < n_ - 1)
                        {
                            grid.connectivity.push_back
                            (
                                localToGlobal[trianglePointIndex(n_, row, col + 1)]
                            );
                            grid.connectivity.push_back
                            (
                                localToGlobal[trianglePointIndex(n_, row + 1, col + 1)]
                            );
                            grid.connectivity.push_back
                            (
                                localToGlobal[trianglePointIndex(n_, row + 1, col)]
                            );
                            grid.offsets.push_back(grid.connectivity.size());
                            grid.cellTypes.push_back(vtkTriType);
                            grid.partitions.push_back(Pstream::myProcNo());
                        }
                    }
                }
            }
            else
            {
                labelList localToGlobal((n_ + 1)*(n_ + 1), -1);

                for (label row = 0; row <= n_; ++row)
                {
                    for (label col = 0; col <= n_; ++col)
                    {
                        const label localPointI = row*(n_ + 1) + col;
                        const vector2D eta
                        (
                            -1.0 + 2.0*scalar(col)/scalar(n_),
                            -1.0 + 2.0*scalar(row)/scalar(n_)
                        );
                        const point physicalPoint =
                            mapToPhysicalQuad(eta, faceVertices);
                        localToGlobal[localPointI] =
                            appendPoint
                            (
                                physicalPoint,
                                ownerCellI,
                                ownerCellVertices,
                                ownerCellType,
                                theta1,
                                theta2
                            );
                    }
                }

                for (label row = 0; row < n_; ++row)
                {
                    for (label col = 0; col < n_; ++col)
                    {
                        const label p0 = row*(n_ + 1) + col;
                        const label p1 = p0 + 1;
                        const label p3 = (row + 1)*(n_ + 1) + col;
                        const label p2 = p3 + 1;

                        grid.connectivity.push_back(localToGlobal[p0]);
                        grid.connectivity.push_back(localToGlobal[p1]);
                        grid.connectivity.push_back(localToGlobal[p2]);
                        grid.connectivity.push_back(localToGlobal[p3]);
                        grid.offsets.push_back(grid.connectivity.size());
                        grid.cellTypes.push_back(vtkQuadType);
                        grid.partitions.push_back(Pstream::myProcNo());
                    }
                }
            }
        }
    }

    const std::vector<dgVtkLagrange::ScalarPointData> scalarArrays
    {
        {"rho", &rhoValues},
        {"E", &EValues},
        {"e", &eValues},
        {"p", &pValues},
        {"T", &TValues},
        {"a", &aValues},
        {"theta1", &theta1Values},
        {"theta2", &theta2Values}
    };

    const std::vector<dgVtkLagrange::VectorPointData> vectorArrays
    {
        {"rhoU", &rhoUValues},
        {"U", &UValues}
    };

    dgVtkLagrange::writeVtuPiece
    (
        filePath,
        runTime_.value(),
        grid,
        scalarArrays,
        vectorArrays
    );
}


void dgFoamSurfaceToVTU::writePvtu(const fileName& filePath) const
{
    const std::vector<dgVtkLagrange::PointDataSchema> pointData
    {
        {"rho", 1},
        {"rhoU", 3},
        {"E", 1},
        {"U", 3},
        {"e", 1},
        {"p", 1},
        {"T", 1},
        {"a", 1},
        {"theta1", 1},
        {"theta2", 1}
    };

    dgVtkLagrange::writePvtu
    (
        filePath,
        runTime_.value(),
        Pstream::nProcs(),
        name(),
        pointData
    );
}

} // End namespace dgFunctionObjects
} // End namespace Foam

// ************************************************************************* //
