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

#include "dgFoamToVTU.H"

#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "IOdictionary.H"
#include "IOobject.H"
#include "DynamicList.H"
#include "cellShape.H"
#include "polyMesh.H"
#include "dgField.H"
#include "dgThermoConservative.H"
#include "dgAffineMapping.H"
#include "boundBox.H"

#include <cmath>
#include <unordered_map>

namespace Foam
{
namespace dgFunctionObjects
{

namespace
{

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


inline label parseAxisName(const word& axisName)
{
    if (axisName == "x")
    {
        return 0;
    }
    if (axisName == "y")
    {
        return 1;
    }
    if (axisName == "z")
    {
        return 2;
    }

    return -1;
}


inline scalar axisAlignment(const vector& dir, const label cartAxis)
{
    return mag(dir[cartAxis])/max(VSMALL, mag(dir));
}


inline vector hexReferenceDirection
(
    const pointField& cellVertices,
    const label refDir
)
{
    switch (refDir)
    {
        case 0:
            return 0.25
               * (
                    (cellVertices[1] - cellVertices[0])
                  + (cellVertices[2] - cellVertices[3])
                  + (cellVertices[5] - cellVertices[4])
                  + (cellVertices[6] - cellVertices[7])
                );

        case 1:
            return 0.25
               * (
                    (cellVertices[3] - cellVertices[0])
                  + (cellVertices[2] - cellVertices[1])
                  + (cellVertices[7] - cellVertices[4])
                  + (cellVertices[6] - cellVertices[5])
                );

        case 2:
            return 0.25
               * (
                    (cellVertices[4] - cellVertices[0])
                  + (cellVertices[5] - cellVertices[1])
                  + (cellVertices[7] - cellVertices[3])
                  + (cellVertices[6] - cellVertices[2])
                );

        default:
            FatalErrorInFunction
                << "Expected HEX reference direction in [0,2], got "
                << refDir << exit(FatalError);
    }

    return vector::zero;
}


inline label collapsedReferenceDirection
(
    const dgCellType type,
    const pointField& cellVertices,
    const label cartAxis
)
{
    switch (type)
    {
        case dgCellType::HEX:
        {
            scalar bestScore = -GREAT;
            label bestDir = -1;

            for (label refDir = 0; refDir < 3; ++refDir)
            {
                const scalar score =
                    axisAlignment
                    (
                        hexReferenceDirection(cellVertices, refDir),
                        cartAxis
                    );

                if (score > bestScore)
                {
                    bestScore = score;
                    bestDir = refDir;
                }
            }

            return bestDir;
        }

        case dgCellType::PRISM:
            return 1;

        default:
            FatalErrorInFunction
                << "Collapsed through-thickness sampling is only implemented "
                << "for HEX and PRISM cells, got " << type
                << exit(FatalError);
    }

    return -1;
}

} // End anonymous namespace

defineTypeNameAndDebug(dgFoamToVTU, 0);
addToRunTimeSelectionTable(functionObject, dgFoamToVTU, dictionary);


dgFoamToVTU::dgFoamToVTU
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
    twoDMode_(false),
    twoDNormalAxis_(-1),
    mergeSharedPoints_(false),
    pointMergeTol_(-1),
    Cp_(Zero),
    Cv_(Zero),
    R_(Zero),
    initialized_(false),
    tessellationPtr_(nullptr)
{
    read(dict);

    Info<< type() << ':' << nl
        << "    " << name << " -> vtu" << nl
        << "        initialized" << endl;
}


bool dgFoamToVTU::read(const dictionary& dict)
{
    functionObject::read(dict);

    // Keep the field names configurable so the function object can be reused
    // with renamed conservative fields or alternate region names.
    regionName_ = dict.getOrDefault<word>("region", polyMesh::defaultRegion);
    rhoName_ = dict.getOrDefault<word>("rho", "rho");
    rhoUName_ = dict.getOrDefault<word>("rhoU", "rhoU");
    EName_ = dict.getOrDefault<word>("E", "E");
    thermoName_ = dict.getOrDefault<word>("thermo", "dgThermoConservative");

    const word twoDNormalAxisName =
        dict.getOrDefault<word>("twoDNormalAxis", word::null);
    twoDMode_ = !twoDNormalAxisName.empty();
    twoDNormalAxis_ = -1;

    if (twoDMode_)
    {
        twoDNormalAxis_ = parseAxisName(twoDNormalAxisName);

        if (twoDNormalAxis_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Expected 'twoDNormalAxis' to be one of x, y, z for "
                << type() << ", got '" << twoDNormalAxisName << "'."
                << exit(FatalIOError);
        }
    }

    mergeSharedPoints_ =
        dict.getOrDefault<Switch>("mergeSharedPoints", false);
    pointMergeTol_ = dict.getOrDefault<scalar>("pointMergeTol", -1);

    if (pointMergeTol_ == 0)
    {
        FatalIOErrorInFunction(dict)
            << "Expected 'pointMergeTol' to be non-zero for dgFoamToVTU."
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
                << "Expected 'outputPoints >= 2' for dgFoamToVTU, got "
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
                << "Expected legacy 'n >= 1' for dgFoamToVTU, got " << n_
                << exit(FatalIOError);
        }
    }
    else
    {
        // Match Nektar's VTU exporter more closely by defaulting to the
        // polynomial's native equi-spaced resolution during initialization.
        n_ = 0;
    }

    pOrder_ = -1;
    Cp_ = Zero;
    Cv_ = Zero;
    R_ = Zero;
    initialized_ = false;
    tessellationPtr_.reset();

    return true;
}


bool dgFoamToVTU::execute()
{
    return true;
}


bool dgFoamToVTU::write()
{
    initialize();

    // Create the standard postProcessing/<name>/<time> directory layout
    // used by OpenFOAM function objects.
    const fileName dir = outputDir();
    mkDir(dir);

    writePiece(pieceFileName());

    if (Pstream::master() && Pstream::parRun())
    {
        writePvtu(pvtuFileName());
    }

    if (log)
    {
        Info<< type() << " wrote VTU samples" << endl;
    }

    return true;
}

void dgFoamToVTU::initialize()
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

    // Read the DG order from the same dgSchemes dictionary used by the solver
    // so the exporter samples with a basis compatible with the active case.
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

    // Pull just the scalar constants needed for pointwise primitive
    // reconstruction. The actual reconstruction logic remains local to this
    // exporter instead of modifying dgThermo itself.
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
            << "Computed Cv = Cp - R <= 0 for dgFoamToVTU."
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

    tessellationPtr_.reset(new dgEquiSpacedTessellation(n_, pOrder_));
    initialized_ = true;
}


const fvMesh& dgFoamToVTU::mesh() const
{
    if (!runTime_.foundObject<fvMesh>(regionName_))
    {
        FatalErrorInFunction
            << "Cannot find fvMesh region '" << regionName_ << "' required by "
            << type() << exit(FatalError);
    }

    return runTime_.lookupObject<fvMesh>(regionName_);
}


const dgField<scalar>& dgFoamToVTU::rhoField() const
{
    return mesh().lookupObject<dgField<scalar>>(rhoName_);
}


const dgField<vector>& dgFoamToVTU::rhoUField() const
{
    return mesh().lookupObject<dgField<vector>>(rhoUName_);
}


const dgField<scalar>& dgFoamToVTU::EField() const
{
    return mesh().lookupObject<dgField<scalar>>(EName_);
}


const dgThermoConservative& dgFoamToVTU::thermo() const
{
    return mesh().lookupObject<dgThermoConservative>(thermoName_);
}


const dgEquiSpacedTessellation& dgFoamToVTU::tessellation() const
{
    if (!tessellationPtr_)
    {
        FatalErrorInFunction
            << type() << " sampling tessellation has not been initialized."
            << exit(FatalError);
    }

    return *tessellationPtr_;
}


fileName dgFoamToVTU::outputDir() const
{
    return
        runTime_.globalPath()
      / functionObject::outputPrefix
      / name()
      / runTime_.timeName();
}


fileName dgFoamToVTU::pieceFileName() const
{
    const fileName base = outputDir();

    if (Pstream::parRun())
    {
        return base/(name() + "_" + Foam::name(Pstream::myProcNo()) + ".vtu");
    }

    return base/(name() + ".vtu");
}


fileName dgFoamToVTU::pvtuFileName() const
{
    return outputDir()/(name() + ".pvtu");
}


template<class Type>
Type dgFoamToVTU::evaluateField
(
    const dgField<Type>& field,
    const label cellI,
    const List<scalar>& basis
)
{
    // Modal DG reconstruction at one arbitrary point:
    //     q(xi) = sum_k q_k * phi_k(xi)
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


void dgFoamToVTU::writePiece(const fileName& filePath) const
{
    const fvMesh& meshRef = mesh();
    const dgField<scalar>& rho = rhoField();
    const dgField<vector>& rhoU = rhoUField();
    const dgField<scalar>& E = EField();

    DynamicList<vector> points;
    DynamicList<scalar> rhoValues;
    DynamicList<vector> rhoUValues;
    DynamicList<scalar> EValues;
    DynamicList<label> sampleCounts;
    DynamicList<label> connectivity;
    DynamicList<label> offsets;
    DynamicList<unsigned char> cellTypes;
    const dgEquiSpacedTessellation& tessellationRef = tessellation();

    points.reserve(meshRef.nCells()*(n_ + 1));
    rhoValues.reserve(meshRef.nCells()*(n_ + 1));
    rhoUValues.reserve(meshRef.nCells()*(n_ + 1));
    EValues.reserve(meshRef.nCells()*(n_ + 1));
    sampleCounts.reserve(meshRef.nCells()*(n_ + 1));

    scalar mergeTol = -1;
    std::unordered_map<PointMergeKey, label, PointMergeKeyHash> mergedPointMap;

    if (mergeSharedPoints_)
    {
        const scalar defaultTol =
            max(SMALL, 1e-12*boundBox(meshRef.points(), false).mag());

        mergeTol = (pointMergeTol_ > 0) ? pointMergeTol_ : defaultTol;
        mergedPointMap.reserve(meshRef.nCells()*(n_ + 1));

        if (log)
        {
            Info<< type()
                << " merging coincident tessellation points with tolerance "
                << mergeTol << endl;
        }
    }

    forAll(meshRef.cells(), cellI)
    {
        const cellShape& shape = meshRef.cellShapes()[cellI];
        const dgEquiSpacedTessellation::CellInfo cellInfo =
            tessellationRef.cellInfo(shape, meshRef.points());
        const pointField& cellVertices = cellInfo.vertices;
        const dgCellType type = cellInfo.type;
        const label collapsedRefDir =
            twoDMode_
          ? collapsedReferenceDirection(type, cellVertices, twoDNormalAxis_)
          : -1;
        const List<dgEquiSpacedTessellation::SamplePoint>& samples =
            twoDMode_
          ? tessellationRef.collapsedSamplePoints(type, collapsedRefDir)
          : tessellationRef.samplePoints(type);
        const label pointOffset = points.size();
        pointField localPoints(samples.size());
        labelList localToGlobal(samples.size(), -1);

        forAll(samples, sampleI)
        {
            const dgEquiSpacedTessellation::SamplePoint& sample = samples[sampleI];

            // Reconstruct conservative variables directly from modal DoFs
            // at the current sample point.
            const scalar rhoValue = evaluateField(rho, cellI, sample.basis);
            const vector rhoUValue = evaluateField(rhoU, cellI, sample.basis);
            const scalar EValue = evaluateField(E, cellI, sample.basis);

            const point physicalPoint = mapEtaToX(sample.eta, cellVertices, type);
            localPoints[sampleI] = physicalPoint;

            if (mergeSharedPoints_)
            {
                const PointMergeKey key = makePointMergeKey(physicalPoint, mergeTol);
                const auto iter = mergedPointMap.find(key);

                if (iter == mergedPointMap.end())
                {
                    const label pointI = points.size();
                    mergedPointMap.emplace(key, pointI);

                    points.append(physicalPoint);
                    rhoValues.append(rhoValue);
                    rhoUValues.append(rhoUValue);
                    EValues.append(EValue);
                    sampleCounts.append(1);

                    localToGlobal[sampleI] = pointI;
                }
                else
                {
                    const label pointI = iter->second;
                    rhoValues[pointI] += rhoValue;
                    rhoUValues[pointI] += rhoUValue;
                    EValues[pointI] += EValue;
                    sampleCounts[pointI] += 1;

                    localToGlobal[sampleI] = pointI;
                }
            }
            else
            {
                points.append(physicalPoint);
                rhoValues.append(rhoValue);
                rhoUValues.append(rhoUValue);
                EValues.append(EValue);
                sampleCounts.append(1);
            }
        }

        if (twoDMode_)
        {
            if (mergeSharedPoints_)
            {
                tessellationRef.appendCollapsedVtkCells
                (
                    type,
                    collapsedRefDir,
                    localPoints,
                    localToGlobal,
                    connectivity,
                    offsets,
                    cellTypes
                );
            }
            else
            {
                tessellationRef.appendCollapsedVtkCells
                (
                    type,
                    collapsedRefDir,
                    localPoints,
                    pointOffset,
                    connectivity,
                    offsets,
                    cellTypes
                );
            }
        }
        else
        {
            if (mergeSharedPoints_)
            {
                tessellationRef.appendVtkCells
                (
                    type,
                    localPoints,
                    localToGlobal,
                    connectivity,
                    offsets,
                    cellTypes
                );
            }
            else
            {
                tessellationRef.appendVtkCells
                (
                    type,
                    localPoints,
                    pointOffset,
                    connectivity,
                    offsets,
                    cellTypes
                );
            }
        }
    }

    List<scalar> pValues(points.size(), Zero);
    List<scalar> TValues(points.size(), Zero);
    List<vector> UValues(points.size(), vector::zero);

    forAll(points, pointI)
    {
        const scalar invCount = 1.0/scalar(sampleCounts[pointI]);
        const scalar rhoValue = rhoValues[pointI]*invCount;
        const vector rhoUValue = rhoUValues[pointI]*invCount;
        const scalar EValue = EValues[pointI]*invCount;

        const scalar rhoSafe =
            (mag(rhoValue) > VSMALL)
          ? rhoValue
          : (rhoValue >= 0.0 ? VSMALL : -VSMALL);

        const vector UValue = rhoUValue/rhoSafe;
        const scalar eValue = EValue/rhoSafe - 0.5*magSqr(UValue);
        const scalar TValue = eValue/Cv_;
        const scalar pValue = rhoValue*R_*TValue;

        pValues[pointI] = pValue;
        TValues[pointI] = TValue;
        UValues[pointI] = UValue;
    }

    mkDir(filePath.path());

    OFstream os(filePath);
    os.precision(16);

    const label nPoints = points.size();
    const label nCells = offsets.size();

    // Write the tessellated sampled solution as a standard VTK XML
    // unstructured grid with linear sub-cells.
    os  << "<?xml version=\"1.0\"?>" << nl
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << nl
        << "  <UnstructuredGrid>" << nl
        << "    <Piece NumberOfPoints=\"" << nPoints
        << "\" NumberOfCells=\"" << nCells << "\">" << nl
        << "      <PointData>" << nl
        << "        <DataArray type=\"Float64\" Name=\"p\" "
        << "NumberOfComponents=\"1\" format=\"ascii\">" << nl;

    forAll(pValues, pointI)
    {
        os << "          " << pValues[pointI] << nl;
    }

    os  << "        </DataArray>" << nl
        << "        <DataArray type=\"Float64\" Name=\"T\" "
        << "NumberOfComponents=\"1\" format=\"ascii\">" << nl;

    forAll(TValues, pointI)
    {
        os << "          " << TValues[pointI] << nl;
    }

    os  << "        </DataArray>" << nl
        << "        <DataArray type=\"Float64\" Name=\"U\" "
        << "NumberOfComponents=\"3\" format=\"ascii\">" << nl;

    forAll(UValues, pointI)
    {
        const vector& value = UValues[pointI];
        os  << "          "
            << value.x() << ' '
            << value.y() << ' '
            << value.z() << nl;
    }

    os  << "        </DataArray>" << nl
        << "      </PointData>" << nl
        << "      <CellData/>" << nl
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


void dgFoamToVTU::writePvtu(const fileName& filePath) const
{
    // The PVTU file is just a lightweight aggregator that points ParaView to
    // all rank-local VTU pieces produced in a parallel run.
    mkDir(filePath.path());

    OFstream os(filePath);

    os  << "<?xml version=\"1.0\"?>" << nl
        << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << nl
        << "  <PUnstructuredGrid GhostLevel=\"0\">" << nl
        << "    <PPointData>" << nl
        << "      <PDataArray type=\"Float64\" Name=\"p\" NumberOfComponents=\"1\"/>" << nl
        << "      <PDataArray type=\"Float64\" Name=\"T\" NumberOfComponents=\"1\"/>" << nl
        << "      <PDataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"3\"/>" << nl
        << "    </PPointData>" << nl
        << "    <PCellData/>" << nl
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
