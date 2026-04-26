/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgFoamVolumeToVTU.H"

#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "polyMesh.H"
#include "dgCellType.H"
#include "dgField.H"
#include "dgAffineMapping.H"
#include "dgGeomMesh.H"
#include "dgThermoConservative.H"
#include "dgVtkLagrangeTools.H"
#include "volFields.H"

#include <vector>

namespace Foam
{
namespace dgFunctionObjects
{

defineTypeNameAndDebug(dgFoamVolumeToVTU, 0);
addToRunTimeSelectionTable(functionObject, dgFoamVolumeToVTU, dictionary);


namespace
{

struct ExportPointData
{
    std::vector<scalar> rho;
    std::vector<vector> rhoU;
    std::vector<scalar> E;
    std::vector<vector> U;
    std::vector<scalar> e;
    std::vector<scalar> p;
    std::vector<scalar> T;
    std::vector<scalar> a;
    std::vector<scalar> theta1;
    std::vector<scalar> theta2;
};


template<class Type>
Type evaluateModalField
(
    const dgField<Type>& field,
    const label cellI,
    const List<scalar>& basis,
    const scalar highOrderScale = 1
)
{
    const cellDof<Type>& cellModes = field.dof()[cellI];

    if (cellModes.nDof() != basis.size())
    {
        FatalErrorInFunction
            << "Basis size mismatch while evaluating field " << field.name()
            << " in cell " << cellI
            << ". DoFs = " << cellModes.nDof()
            << ", basis size = " << basis.size()
            << exit(FatalError);
    }

    Type value = pTraits<Type>::zero;

    for (label modeI = 0; modeI < cellModes.nDof(); ++modeI)
    {
        const scalar modeScale = (modeI == 0 ? scalar(1) : highOrderScale);
        value += modeScale*basis[modeI]*cellModes[modeI];
    }

    return value;
}

} // End anonymous namespace


dgFoamVolumeToVTU::dgFoamVolumeToVTU
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
    initialized_(false)
{
    read(dict);

    Info<< type() << ':' << nl
        << "    " << name << " -> native VTK Lagrange vtu" << nl
        << "        initialized" << endl;
}


bool dgFoamVolumeToVTU::read(const dictionary& dict)
{
    functionObject::read(dict);

    regionName_ = dict.getOrDefault<word>("region", polyMesh::defaultRegion);
    rhoName_ = dict.getOrDefault<word>("rho", "rho");
    rhoUName_ = dict.getOrDefault<word>("rhoU", "rhoU");
    EName_ = dict.getOrDefault<word>("E", "E");
    thermoName_ = dict.getOrDefault<word>("thermo", "dgThermoConservative");

    const wordList ignoredKeys
    (
        {
            "outputPoints",
            "n",
            "twoDNormalAxis",
            "mergeSharedPoints",
            "pointMergeTol"
        }
    );

    for (const word& key : ignoredKeys)
    {
        if (dict.found(key))
        {
            WarningInFunction
                << type() << " now writes native VTK Lagrange cells and "
                << "ignores legacy option '" << key << "'." << endl;
        }
    }

    initialized_ = false;
    return true;
}


bool dgFoamVolumeToVTU::execute()
{
    return true;
}


bool dgFoamVolumeToVTU::write()
{
    initialize();
    mkDir(outputDir());

    writePiece(pieceFileName());

    if (Pstream::master() && Pstream::parRun())
    {
        writePvtu(pvtuFileName());
    }

    if (log)
    {
        Info<< type() << " wrote native VTK Lagrange samples" << endl;
    }

    return true;
}


void dgFoamVolumeToVTU::initialize()
{
    if (initialized_)
    {
        return;
    }

    const dgThermoConservative& dgThermo = thermo();

    if (!dgThermo.heIsInternalEnergy())
    {
        FatalErrorInFunction
            << type()
            << " reconstructs primitive variables from conservative energy "
            << "assuming internal-energy storage." << nl
            << "Detected thermo energy variable: " << dgThermo.heName()
            << exit(FatalError);
    }

    initialized_ = true;
}


const fvMesh& dgFoamVolumeToVTU::mesh() const
{
    if (!runTime_.foundObject<fvMesh>(regionName_))
    {
        FatalErrorInFunction
            << "Cannot find fvMesh region '" << regionName_ << "' required by "
            << type() << exit(FatalError);
    }

    return runTime_.lookupObject<fvMesh>(regionName_);
}


const dgField<scalar>& dgFoamVolumeToVTU::rhoField() const
{
    return mesh().lookupObject<dgField<scalar>>(rhoName_);
}


const dgField<vector>& dgFoamVolumeToVTU::rhoUField() const
{
    return mesh().lookupObject<dgField<vector>>(rhoUName_);
}


const dgField<scalar>& dgFoamVolumeToVTU::EField() const
{
    return mesh().lookupObject<dgField<scalar>>(EName_);
}


const dgThermoConservative& dgFoamVolumeToVTU::thermo() const
{
    return mesh().lookupObject<dgThermoConservative>(thermoName_);
}


fileName dgFoamVolumeToVTU::outputDir() const
{
    return
        runTime_.globalPath()
      / functionObject::outputPrefix
      / name()
      / runTime_.timeName();
}


fileName dgFoamVolumeToVTU::pieceFileName() const
{
    const fileName base = outputDir();

    if (Pstream::parRun())
    {
        return base/(name() + "_" + Foam::name(Pstream::myProcNo()) + ".vtu");
    }

    return base/(name() + ".vtu");
}


fileName dgFoamVolumeToVTU::pvtuFileName() const
{
    return outputDir()/(name() + ".pvtu");
}


void dgFoamVolumeToVTU::writePiece(const fileName& filePath) const
{
    const fvMesh& meshRef = mesh();
    const dgField<scalar>& rhoFieldRef = rhoField();
    const dgField<vector>& rhoUFieldRef = rhoUField();
    const dgField<scalar>& EFieldRef = EField();
    const dgThermoConservative& thermoRef = thermo();

    const volScalarField* theta1FieldPtr =
        meshRef.foundObject<volScalarField>("theta1")
      ? &meshRef.lookupObject<volScalarField>("theta1")
      : nullptr;

    const volScalarField* theta2FieldPtr =
        meshRef.foundObject<volScalarField>("theta2")
      ? &meshRef.lookupObject<volScalarField>("theta2")
      : nullptr;

    dgGeomMesh dgMesh(meshRef);

    const label pOrder = dgMesh.pOrder();
    const label hoOrder = dgVtkLagrange::exportOrder(pOrder);

    if (hoOrder > 8)
    {
        FatalErrorInFunction
            << "Current native VTK Lagrange exporter supports order up to 8. "
            << "Detected DG pOrder = " << pOrder
            << ", which maps to VTK order " << hoOrder
            << exit(FatalError);
    }

    dgVtkLagrange::GridBuffers grid;
    ExportPointData pointData;
    label skippedPyramids = 0;

    for (label cellI = 0; cellI < dgMesh.nCells(); ++cellI)
    {
        const dgGeomCell* cellPtr = dgMesh.cells()[cellI];

        if (!cellPtr)
        {
            continue;
        }

        const dgCellType type = cellPtr->type();

        if (!dgVtkLagrange::supportsNativeLagrange(type))
        {
            if (type == dgCellType::PYRAMID)
            {
                ++skippedPyramids;
                continue;
            }

            FatalErrorInFunction
                << "Unsupported DG cell type "
                << dgVtkLagrange::cellTypeName(type)
                << " in cell " << cellI << exit(FatalError);
        }

        const pointField& cellVertices = cellPtr->points();
        const std::vector<vector> pyfrPts =
            dgVtkLagrange::pyfrStdElementPoints(type, hoOrder);
        const std::vector<label> nodemap =
            dgVtkLagrange::vtkNodemap(type, pyfrPts.size());
        const label pointOffset = grid.points.size();
        const scalar theta1 =
            theta1FieldPtr ? theta1FieldPtr->primitiveField()[cellI] : scalar(1);
        const scalar theta2 =
            theta2FieldPtr ? theta2FieldPtr->primitiveField()[cellI] : scalar(1);
        const scalar rhoTheta = theta1*theta2;

        for (const label pyfrPointI : nodemap)
        {
            const vector xi = pyfrPts[pyfrPointI];
            const vector eta = dgVtkLagrange::xiToEta(type, xi);

            List<scalar> basis;
            dgVtkLagrange::computeBasisAt(type, eta, pOrder, basis);

            const point x = mapEtaToX(eta, cellVertices, type);
            const scalar rhoVal =
                evaluateModalField(rhoFieldRef, cellI, basis, rhoTheta);
            const vector rhoUVal =
                evaluateModalField(rhoUFieldRef, cellI, basis, theta2);
            const scalar EVal =
                evaluateModalField(EFieldRef, cellI, basis, theta2);

            const scalar rhoSafe = max(rhoVal, SMALL);
            const vector UVal = rhoUVal/rhoSafe;
            const scalar eVal = EVal/rhoSafe - 0.5*magSqr(UVal);

            grid.points.push_back(x);
            pointData.rho.push_back(rhoVal);
            pointData.rhoU.push_back(rhoUVal);
            pointData.E.push_back(EVal);
            pointData.U.push_back(UVal);
            pointData.e.push_back(eVal);
            pointData.T.push_back(thermoRef.calcTemperatureFromRhoHe(cellI, rhoSafe, eVal));
            pointData.p.push_back(thermoRef.calcPressureFromRhoHe(cellI, rhoSafe, eVal));
            pointData.a.push_back(thermoRef.calcSpeedOfSoundFromRhoHe(cellI, rhoSafe, eVal));
            pointData.theta1.push_back(theta1);
            pointData.theta2.push_back(theta2);
        }

        for (label localPointI = 0; localPointI < label(nodemap.size()); ++localPointI)
        {
            grid.connectivity.push_back(pointOffset + localPointI);
        }

        grid.offsets.push_back(grid.connectivity.size());
        grid.cellTypes.push_back(dgVtkLagrange::vtkLagrangeType(type));
        grid.partitions.push_back(Pstream::myProcNo());
    }

    const std::vector<dgVtkLagrange::ScalarPointData> scalarArrays
    {
        {"rho", &pointData.rho},
        {"E", &pointData.E},
        {"e", &pointData.e},
        {"p", &pointData.p},
        {"T", &pointData.T},
        {"a", &pointData.a},
        {"theta1", &pointData.theta1},
        {"theta2", &pointData.theta2}
    };

    const std::vector<dgVtkLagrange::VectorPointData> vectorArrays
    {
        {"rhoU", &pointData.rhoU},
        {"U", &pointData.U}
    };

    dgVtkLagrange::writeVtuPiece
    (
        filePath,
        runTime_.value(),
        grid,
        scalarArrays,
        vectorArrays
    );

    if (skippedPyramids > 0)
    {
        WarningInFunction
            << "Skipped " << skippedPyramids
            << " pyramid cells because native VTK_LAGRANGE pyramid export is "
            << "not implemented yet." << endl;
    }
}


void dgFoamVolumeToVTU::writePvtu(const fileName& filePath) const
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
