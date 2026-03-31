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

namespace Foam
{
namespace dgFunctionObjects
{

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
    n_(1),
    pOrder_(-1),
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

    n_ = readLabel(dict.lookup("n"));

    if (n_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Expected 'n >= 1' for dgFoamToVTU, got " << n_
            << exit(FatalIOError);
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
    DynamicList<scalar> pValues;
    DynamicList<scalar> TValues;
    DynamicList<vector> UValues;
    DynamicList<label> connectivity;
    DynamicList<label> offsets;
    DynamicList<unsigned char> cellTypes;
    const dgEquiSpacedTessellation& tessellationRef = tessellation();

    points.reserve(meshRef.nCells()*(n_ + 1));
    pValues.reserve(meshRef.nCells()*(n_ + 1));
    TValues.reserve(meshRef.nCells()*(n_ + 1));
    UValues.reserve(meshRef.nCells()*(n_ + 1));

    forAll(meshRef.cells(), cellI)
    {
        const cellShape& shape = meshRef.cellShapes()[cellI];
        const dgEquiSpacedTessellation::CellInfo cellInfo =
            tessellationRef.cellInfo(shape, meshRef.points());
        const pointField& cellVertices = cellInfo.vertices;
        const dgCellType type = cellInfo.type;
        const List<dgEquiSpacedTessellation::SamplePoint>& samples =
            tessellationRef.samplePoints(type);
        const label pointOffset = points.size();
        pointField localPoints(samples.size());

        forAll(samples, sampleI)
        {
            const dgEquiSpacedTessellation::SamplePoint& sample = samples[sampleI];

            // Reconstruct conservative variables directly from modal DoFs
            // at the current sample point.
            const scalar rhoValue = evaluateField(rho, cellI, sample.basis);
            const vector rhoUValue = evaluateField(rhoU, cellI, sample.basis);
            const scalar EValue = evaluateField(E, cellI, sample.basis);

            const scalar rhoSafe =
                (mag(rhoValue) > VSMALL)
              ? rhoValue
              : (rhoValue >= 0.0 ? VSMALL : -VSMALL);

            const vector UValue = rhoUValue/rhoSafe;
            const scalar eValue = EValue/rhoSafe - 0.5*magSqr(UValue);
            const scalar TValue = eValue/Cv_;
            const scalar pValue = rhoValue*R_*TValue;

            // Store the physical-space point and the derived primitive data
            // as point data in the VTU file.
            const point physicalPoint = mapEtaToX(sample.eta, cellVertices, type);

            localPoints[sampleI] = physicalPoint;
            points.append(physicalPoint);
            pValues.append(pValue);
            TValues.append(TValue);
            UValues.append(UValue);
        }

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
