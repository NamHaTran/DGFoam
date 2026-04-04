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
#include "basisFunctions.H"
#include "boundBox.H"
#include "refCoordTransforms.H"

#include <cmath>
#include <initializer_list>
#include <unordered_map>

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


inline void appendVtkCell
(
    std::initializer_list<label> pointLabels,
    const unsigned char vtkType,
    DynamicList<label>& connectivity,
    DynamicList<label>& offsets,
    DynamicList<unsigned char>& cellTypes
)
{
    for (const label pointI : pointLabels)
    {
        connectivity.append(pointI);
    }

    offsets.append(connectivity.size());
    cellTypes.append(vtkType);
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
    Cp_(Zero),
    Cv_(Zero),
    R_(Zero),
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
    Cp_ = Zero;
    Cv_ = Zero;
    R_ = Zero;
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

    if
    (
        dgThermo.eos().type() != "idealGas"
     || dgThermo.thermo().type() != "constantCp"
     || dgThermo.energyModel().type() != "sensibleInternalEnergy"
    )
    {
        FatalErrorInFunction
            << type() << " currently supports only the thermo chain" << nl
            << "  eos    = idealGas" << nl
            << "  thermo = constantCp" << nl
            << "  energy = sensibleInternalEnergy" << nl
            << "Detected: eos=" << dgThermo.eos().type()
            << ", thermo=" << dgThermo.thermo().type()
            << ", energy=" << dgThermo.energyModel().type()
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

    IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime_.constant(),
            meshRef,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (!thermoDict.found("mixture"))
    {
        FatalIOErrorInFunction(thermoDict)
            << "Missing 'mixture' sub-dictionary in thermophysicalProperties."
            << exit(FatalIOError);
    }

    const dictionary& mixDict = thermoDict.subDict("mixture");

    if (!mixDict.found("thermodynamics"))
    {
        FatalIOErrorInFunction(mixDict)
            << "Missing 'thermodynamics' sub-dictionary in 'mixture' for "
            << type() << '.'
            << exit(FatalIOError);
    }

    const dictionary& thermoCoeffs = mixDict.subDict("thermodynamics");
    Cp_ = readScalar(thermoCoeffs.lookup("Cp"));
    R_ = dgThermo.R();
    Cv_ = Cp_ - R_;

    if (Cv_ <= SMALL)
    {
        FatalIOErrorInFunction(thermoDict)
            << "Computed Cv = Cp - R <= 0 for dgFoamSurfaceToVTU."
            << nl
            << "Cp = " << Cp_ << ", R = " << R_
            << exit(FatalIOError);
    }

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
    const List<scalar>& basis
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
        value += cellDof[dofI]*basis[dofI];
    }

    return value;
}


void dgFoamSurfaceToVTU::writePiece(const fileName& filePath) const
{
    const fvMesh& meshRef = mesh();
    const dgField<scalar>& rho = rhoField();
    const dgField<vector>& rhoU = rhoUField();
    const dgField<scalar>& E = EField();

    DynamicList<vector> points;
    DynamicList<scalar> rhoCellValues;
    DynamicList<scalar> pCellValues;
    DynamicList<scalar> TCellValues;
    DynamicList<vector> UCellValues;
    DynamicList<label> connectivity;
    DynamicList<label> offsets;
    DynamicList<unsigned char> cellTypes;

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

            if (nFacePoints == 3)
            {
                labelList localToGlobal(trianglePointCount(n_), -1);
                pointField localPoints(trianglePointCount(n_), point::zero);

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
                        localPoints[localPointI] = physicalPoint;

                        if (mergeSharedPoints_)
                        {
                            const PointMergeKey key =
                                makePointMergeKey(physicalPoint, mergeTol);
                            const auto iter = mergedPointMap.find(key);

                            if (iter == mergedPointMap.end())
                            {
                                const label pointI = points.size();
                                mergedPointMap.emplace(key, pointI);

                                points.append(physicalPoint);
                                localToGlobal[localPointI] = pointI;
                            }
                            else
                            {
                                const label pointI = iter->second;
                                localToGlobal[localPointI] = pointI;
                            }
                        }
                        else
                        {
                            const label pointI = points.size();
                            points.append(physicalPoint);
                            localToGlobal[localPointI] = pointI;
                        }
                    }
                }

                for (label row = 0; row < n_; ++row)
                {
                    for (label col = 0; col < n_ - row; ++col)
                    {
                        appendVtkCell
                        (
                            {
                                localToGlobal[trianglePointIndex(n_, row, col)],
                                localToGlobal[trianglePointIndex(n_, row, col + 1)],
                                localToGlobal[trianglePointIndex(n_, row + 1, col)]
                            },
                            vtkTriType,
                            connectivity,
                            offsets,
                            cellTypes
                        );

                        {
                            const label p0 = trianglePointIndex(n_, row, col);
                            const label p1 = trianglePointIndex(n_, row, col + 1);
                            const label p2 = trianglePointIndex(n_, row + 1, col);
                            const point subFaceCenter =
                                (localPoints[p0] + localPoints[p1] + localPoints[p2])/3.0;
                            const vector ownerEta =
                                mapXToEta(subFaceCenter, ownerCellVertices, ownerCellType);
                            List<scalar> basis;
                            computeBasisAt(ownerCellType, ownerEta, pOrder_, basis);

                            const scalar rhoValue =
                                evaluateField(rho, ownerCellI, basis);
                            const vector rhoUValue =
                                evaluateField(rhoU, ownerCellI, basis);
                            const scalar EValue =
                                evaluateField(E, ownerCellI, basis);
                            const scalar rhoSafe =
                                (mag(rhoValue) > VSMALL)
                              ? rhoValue
                              : (rhoValue >= 0.0 ? VSMALL : -VSMALL);
                            const vector UValue = rhoUValue/rhoSafe;
                            const scalar eValue =
                                EValue/rhoSafe - 0.5*magSqr(UValue);
                            const scalar TValue = eValue/Cv_;
                            const scalar pValue = rhoValue*R_*TValue;

                            rhoCellValues.append(rhoValue);
                            pCellValues.append(pValue);
                            TCellValues.append(TValue);
                            UCellValues.append(UValue);
                        }

                        if (row + col < n_ - 1)
                        {
                            appendVtkCell
                            (
                                {
                                    localToGlobal[trianglePointIndex(n_, row, col + 1)],
                                    localToGlobal[trianglePointIndex(n_, row + 1, col + 1)],
                                    localToGlobal[trianglePointIndex(n_, row + 1, col)]
                                },
                                vtkTriType,
                                connectivity,
                                offsets,
                                cellTypes
                            );

                            {
                                const label p0 = trianglePointIndex(n_, row, col + 1);
                                const label p1 = trianglePointIndex(n_, row + 1, col + 1);
                                const label p2 = trianglePointIndex(n_, row + 1, col);
                                const point subFaceCenter =
                                    (localPoints[p0] + localPoints[p1] + localPoints[p2])/3.0;
                                const vector ownerEta =
                                    mapXToEta(subFaceCenter, ownerCellVertices, ownerCellType);
                                List<scalar> basis;
                                computeBasisAt(ownerCellType, ownerEta, pOrder_, basis);

                                const scalar rhoValue =
                                    evaluateField(rho, ownerCellI, basis);
                                const vector rhoUValue =
                                    evaluateField(rhoU, ownerCellI, basis);
                                const scalar EValue =
                                    evaluateField(E, ownerCellI, basis);
                                const scalar rhoSafe =
                                    (mag(rhoValue) > VSMALL)
                                  ? rhoValue
                                  : (rhoValue >= 0.0 ? VSMALL : -VSMALL);
                                const vector UValue = rhoUValue/rhoSafe;
                                const scalar eValue =
                                    EValue/rhoSafe - 0.5*magSqr(UValue);
                                const scalar TValue = eValue/Cv_;
                                const scalar pValue = rhoValue*R_*TValue;

                                rhoCellValues.append(rhoValue);
                                pCellValues.append(pValue);
                                TCellValues.append(TValue);
                                UCellValues.append(UValue);
                            }
                        }
                    }
                }
            }
            else
            {
                labelList localToGlobal((n_ + 1)*(n_ + 1), -1);
                pointField localPoints((n_ + 1)*(n_ + 1), point::zero);

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
                        localPoints[localPointI] = physicalPoint;

                        if (mergeSharedPoints_)
                        {
                            const PointMergeKey key =
                                makePointMergeKey(physicalPoint, mergeTol);
                            const auto iter = mergedPointMap.find(key);

                            if (iter == mergedPointMap.end())
                            {
                                const label pointI = points.size();
                                mergedPointMap.emplace(key, pointI);

                                points.append(physicalPoint);
                                localToGlobal[localPointI] = pointI;
                            }
                            else
                            {
                                const label pointI = iter->second;
                                localToGlobal[localPointI] = pointI;
                            }
                        }
                        else
                        {
                            const label pointI = points.size();
                            points.append(physicalPoint);
                            localToGlobal[localPointI] = pointI;
                        }
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

                        appendVtkCell
                        (
                            {
                                localToGlobal[p0],
                                localToGlobal[p1],
                                localToGlobal[p2],
                                localToGlobal[p3]
                            },
                            vtkQuadType,
                            connectivity,
                            offsets,
                            cellTypes
                        );

                        {
                            const point subFaceCenter =
                                (localPoints[p0] + localPoints[p1]
                               + localPoints[p2] + localPoints[p3])/4.0;
                            const vector ownerEta =
                                mapXToEta(subFaceCenter, ownerCellVertices, ownerCellType);
                            List<scalar> basis;
                            computeBasisAt(ownerCellType, ownerEta, pOrder_, basis);

                            const scalar rhoValue =
                                evaluateField(rho, ownerCellI, basis);
                            const vector rhoUValue =
                                evaluateField(rhoU, ownerCellI, basis);
                            const scalar EValue =
                                evaluateField(E, ownerCellI, basis);
                            const scalar rhoSafe =
                                (mag(rhoValue) > VSMALL)
                              ? rhoValue
                              : (rhoValue >= 0.0 ? VSMALL : -VSMALL);
                            const vector UValue = rhoUValue/rhoSafe;
                            const scalar eValue =
                                EValue/rhoSafe - 0.5*magSqr(UValue);
                            const scalar TValue = eValue/Cv_;
                            const scalar pValue = rhoValue*R_*TValue;

                            rhoCellValues.append(rhoValue);
                            pCellValues.append(pValue);
                            TCellValues.append(TValue);
                            UCellValues.append(UValue);
                        }
                    }
                }
            }
        }
    }

    mkDir(filePath.path());

    OFstream os(filePath);
    os.precision(16);

    const label nPoints = points.size();
    const label nCells = offsets.size();

    if
    (
        rhoCellValues.size() != nCells
     || pCellValues.size() != nCells
     || TCellValues.size() != nCells
     || UCellValues.size() != nCells
    )
    {
        FatalErrorInFunction
            << "Cell-data count mismatch while writing surface VTU: nCells="
            << nCells
            << ", rho=" << rhoCellValues.size()
            << ", p=" << pCellValues.size()
            << ", T=" << TCellValues.size()
            << ", U=" << UCellValues.size()
            << exit(FatalError);
    }

    os  << "<?xml version=\"1.0\"?>" << nl
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << nl
        << "  <UnstructuredGrid>" << nl
        << "    <Piece NumberOfPoints=\"" << nPoints
        << "\" NumberOfCells=\"" << nCells << "\">" << nl
        << "      <PointData/>" << nl
        << "      <CellData>" << nl
        << "        <DataArray type=\"Float64\" Name=\"rho\" "
        << "NumberOfComponents=\"1\" format=\"ascii\">" << nl;

    forAll(rhoCellValues, cellI)
    {
        os << "          " << rhoCellValues[cellI] << nl;
    }

    os  << "        </DataArray>" << nl
        << "        <DataArray type=\"Float64\" Name=\"p\" "
        << "NumberOfComponents=\"1\" format=\"ascii\">" << nl;

    forAll(pCellValues, cellI)
    {
        os << "          " << pCellValues[cellI] << nl;
    }

    os  << "        </DataArray>" << nl
        << "        <DataArray type=\"Float64\" Name=\"T\" "
        << "NumberOfComponents=\"1\" format=\"ascii\">" << nl;

    forAll(TCellValues, cellI)
    {
        os << "          " << TCellValues[cellI] << nl;
    }

    os  << "        </DataArray>" << nl
        << "        <DataArray type=\"Float64\" Name=\"U\" "
        << "NumberOfComponents=\"3\" format=\"ascii\">" << nl;

    forAll(UCellValues, cellI)
    {
        const vector& value = UCellValues[cellI];
        os  << "          "
            << value.x() << ' '
            << value.y() << ' '
            << value.z() << nl;
    }

    os  << "        </DataArray>" << nl
        << "      </CellData>" << nl
        << "      <Points>" << nl
        << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        << "format=\"ascii\">" << nl;

    forAll(points, pointI)
    {
        const vector& point = points[pointI];
        os  << "          "
            << point.x() << ' '
            << point.y() << ' '
            << point.z() << nl;
    }

    os  << "        </DataArray>" << nl
        << "      </Points>" << nl
        << "      <Cells>" << nl
        << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
        << "format=\"ascii\">" << nl;

    forAll(connectivity, connI)
    {
        os << "          " << connectivity[connI] << nl;
    }

    os  << "        </DataArray>" << nl
        << "        <DataArray type=\"Int64\" Name=\"offsets\" "
        << "format=\"ascii\">" << nl;

    forAll(offsets, cellI)
    {
        os << "          " << offsets[cellI] << nl;
    }

    os  << "        </DataArray>" << nl
        << "        <DataArray type=\"UInt8\" Name=\"types\" "
        << "format=\"ascii\">" << nl;

    forAll(cellTypes, cellI)
    {
        os << "          " << label(cellTypes[cellI]) << nl;
    }

    os  << "        </DataArray>" << nl
        << "      </Cells>" << nl
        << "    </Piece>" << nl
        << "  </UnstructuredGrid>" << nl
        << "</VTKFile>" << nl;
}


void dgFoamSurfaceToVTU::writePvtu(const fileName& filePath) const
{
    mkDir(filePath.path());

    OFstream os(filePath);

    os  << "<?xml version=\"1.0\"?>" << nl
        << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << nl
        << "  <PUnstructuredGrid GhostLevel=\"0\">" << nl
        << "    <PPointData/>" << nl
        << "    <PCellData>" << nl
        << "      <PDataArray type=\"Float64\" Name=\"rho\" NumberOfComponents=\"1\"/>" << nl
        << "      <PDataArray type=\"Float64\" Name=\"p\" NumberOfComponents=\"1\"/>" << nl
        << "      <PDataArray type=\"Float64\" Name=\"T\" NumberOfComponents=\"1\"/>" << nl
        << "      <PDataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"3\"/>" << nl
        << "    </PCellData>" << nl
        << "    <PPoints>" << nl
        << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>" << nl
        << "    </PPoints>" << nl;

    for (label procI = 0; procI < Pstream::nProcs(); ++procI)
    {
        os  << "    <Piece Source=\""
            << name() << "_" << procI << ".vtu\"/>" << nl;
    }

    os  << "  </PUnstructuredGrid>" << nl
        << "</VTKFile>" << nl;
}

} // End namespace dgFunctionObjects
} // End namespace Foam

// ************************************************************************* //
