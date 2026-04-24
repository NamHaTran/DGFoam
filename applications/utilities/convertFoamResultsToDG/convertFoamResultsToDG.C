/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Convert finite-volume compressible fields p/T/U/rho at a selected time to
    DGFoam conservative restart fields rho/rhoU/E. The DG polynomial order is
    read automatically from system/dgSchemes through dgGeomMesh. For pOrder=0
    the restart is piecewise constant; for pOrder>=1 the utility reconstructs
    a linear FV state in each cell and L2-projects it onto the available DG
    basis. Modes above linear order therefore receive the consistent
    projection of a linear polynomial, which is zero for orthogonal bases up to
    quadrature tolerance.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"

#include "Jacobian.H"
#include "dgAffineMapping.H"
#include "dgField.H"
#include "dgGeomMesh.H"
#include "eqnOfState.H"
#include "thermoLaw.H"
#include "energy.H"

#include <fstream>
#include <regex>
#include <sstream>
#include <string>

namespace Foam
{
namespace
{

scalar calcStoredThermoEnergy
(
    const scalar rho,
    const scalar T,
    const eqnOfState& eos,
    const energy& energyModel,
    const bool heIsInternalEnergy
)
{
    if (heIsInternalEnergy && eos.canCalcEFromRhoT())
    {
        return eos.calcEFromRhoT(rho, T);
    }

    return energyModel.calcHe(T);
}

template<class Type>
void setConstantCellDof
(
    dgField<Type>& field,
    const label cellI,
    const Type& value
)
{
    cellDof<Type>& dofs = field.dof()[cellI];

    if (dofs.nDof() == 0)
    {
        FatalErrorInFunction
            << "Field " << field.name() << " has no DoFs in cell "
            << cellI << exit(FatalError);
    }

    dofs[0] = value;

    for (label dofI = 1; dofI < dofs.nDof(); ++dofI)
    {
        dofs[dofI] = pTraits<Type>::zero;
    }

    field.dof().updateCellDof(cellI);
    field.gaussFields()[cellI].interpolateFromDof();
}


void initialiseConservativeSeed
(
    dgField<scalar>& rho,
    dgField<vector>& rhoU,
    dgField<scalar>& E,
    const volScalarField& rhoFoam,
    const volScalarField& EFoam,
    const volVectorField& rhoUFoam
)
{
    forAll(rhoFoam, cellI)
    {
        const scalar rhoValue = rhoFoam[cellI];

        setConstantCellDof(rho, cellI, rhoValue);
        setConstantCellDof(rhoU, cellI, rhoUFoam[cellI]);
        setConstantCellDof(E, cellI, EFoam[cellI]);
    }
}


bool replaceAll
(
    std::string& text,
    const std::string& from,
    const std::string& to
)
{
    bool changed = false;
    std::string::size_type pos = 0;

    while ((pos = text.find(from, pos)) != std::string::npos)
    {
        text.replace(pos, from.size(), to);
        pos += to.size();
        changed = true;
    }

    return changed;
}


bool convertSymmetryPlanePatchType(const fileName& filePath)
{
    std::ifstream input(filePath.c_str());

    if (!input.good())
    {
        return false;
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();

    std::string content = buffer.str();
    input.close();

    const bool convertedSymmetryPlane =
        replaceAll(content, "symmetryPlane", "symmetry");

    const std::regex emptyPatchType
    (
        R"((\btype[ \t]+)empty([ \t]*;))"
    );
    const std::string convertedContent =
        std::regex_replace(content, emptyPatchType, "$1symmetry$2");

    const bool convertedEmpty = (convertedContent != content);
    content = convertedContent;

    if (!convertedSymmetryPlane && !convertedEmpty)
    {
        return false;
    }

    std::ofstream output(filePath.c_str(), std::ios::trunc);
    output << content;

    if (!output.good())
    {
        FatalErrorInFunction
            << "Failed to rewrite " << filePath
            << " after converting symmetryPlane to symmetry."
            << exit(FatalError);
    }

    return true;
}


label maxDofPerCell(const dgGeomMesh& dgMesh)
{
    label maxDof = 0;

    for (label cellI = 0; cellI < dgMesh.nCells(); ++cellI)
    {
        maxDof = max(maxDof, label(dgMesh.cells()[cellI]->nDof()));
    }

    return maxDof;
}


void convertWrittenSymmetryPlanePatchTypes
(
    const Time& runTime,
    const dgGeomMesh& dgMesh
)
{
    const label maxDof = maxDofPerCell(dgMesh);
    label nConverted = 0;

    wordList fieldNames(3);
    fieldNames[0] = "rho";
    fieldNames[1] = "rhoU";
    fieldNames[2] = "E";

    forAll(fieldNames, fieldI)
    {
        const word& fieldName = fieldNames[fieldI];

        if (convertSymmetryPlanePatchType(runTime.timePath()/fieldName))
        {
            ++nConverted;
        }

        for (label dofI = 0; dofI < maxDof; ++dofI)
        {
            const word dofFieldName =
                fieldName + "_dof" + Foam::name(dofI);

            if
            (
                convertSymmetryPlanePatchType
                (
                    runTime.timePath()/dofFieldName
                )
            )
            {
                ++nConverted;
            }
        }
    }

    if (nConverted)
    {
        Info<< "Converted symmetryPlane/empty patch-field entries to symmetry in "
            << nConverted << " written DG restart files." << nl;
    }
}


void convertInputFoamPatchTypes(const Time& runTime)
{
    label nConverted = 0;

    wordList fieldNames(4);
    fieldNames[0] = "p";
    fieldNames[1] = "T";
    fieldNames[2] = "U";
    fieldNames[3] = "rho";

    forAll(fieldNames, fieldI)
    {
        if (convertSymmetryPlanePatchType(runTime.timePath()/fieldNames[fieldI]))
        {
            ++nConverted;
        }
    }

    if (nConverted)
    {
        Info<< "Converted symmetryPlane/empty patch-field entries to symmetry in "
            << nConverted << " input finite-volume files." << nl;
    }
}


void buildBoundaryConservativeFields
(
    volScalarField& rhoFoam,
    volVectorField& rhoUFoam,
    volScalarField& EFoam,
    const volScalarField& pFoam,
    const volScalarField& TFoam,
    const volVectorField& UFoam,
    const eqnOfState& eos,
    const energy& energyModel,
    const bool heIsInternalEnergy
)
{
    forAll(rhoUFoam.boundaryField(), patchI)
    {
        auto& rhoPatch = rhoFoam.boundaryFieldRef()[patchI];
        auto& rhoUPatch = rhoUFoam.boundaryFieldRef()[patchI];
        auto& EPatch = EFoam.boundaryFieldRef()[patchI];

        const auto& pPatch = pFoam.boundaryField()[patchI];
        const auto& TPatch = TFoam.boundaryField()[patchI];
        const auto& UPatch = UFoam.boundaryField()[patchI];
        const labelUList& faceCells = rhoFoam.boundaryField()[patchI].patch().faceCells();

        forAll(faceCells, faceI)
        {
            const label cellI = faceCells[faceI];
            const scalar rhoValue = eos.calcRhoFromPT(pPatch[faceI], TPatch[faceI]);
            const scalar TValue = TPatch[faceI];
            const vector UValue = UPatch[faceI];
            const scalar heValue =
                calcStoredThermoEnergy
                (
                    rhoValue,
                    TValue,
                    eos,
                    energyModel,
                    heIsInternalEnergy
                );

            rhoPatch[faceI] = rhoValue;
            rhoUPatch[faceI] = rhoValue*UValue;
            EPatch[faceI] = rhoValue*(heValue + 0.5*magSqr(UValue));
        }
    }
}


template<class Type, class GradType>
Type evaluateLinearReconstruction
(
    const Type& meanValue,
    const GradType& gradValue,
    const vector& dX
)
{
    return meanValue + (gradValue & dX);
}


template<class Type, class GradType>
void projectLinearCellField
(
    dgField<Type>& field,
    const label cellI,
    const UList<Type>& meanField,
    const UList<GradType>& gradField,
    const dgGeomMesh& dgMesh
)
{
    const dgGeomCell& cell = *dgMesh.cells()[cellI];
    const List<vector>& gaussEta = cell.gaussPoints();
    const List<scalar>& weights = cell.weights();
    const List<List<scalar>>& basis = cell.basis();
    const List<scalar> massDiag = cell.massMatrixDiag();
    const List<vector>& cellPoints = cell.points();
    const point xC = cell.centre();

    cellDof<Type>& cellModes = field.dof()[cellI];

    if (cellModes.nDof() != massDiag.size())
    {
        FatalErrorInFunction
            << "DoF count mismatch for field " << field.name()
            << " in cell " << cellI << ": dof=" << cellModes.nDof()
            << ", massDiag=" << massDiag.size() << exit(FatalError);
    }

    List<Type> rhs(cellModes.nDof(), pTraits<Type>::zero);

    forAll(gaussEta, gpI)
    {
        const vector xGp =
            Foam::mapEtaToX(gaussEta[gpI], cellPoints, cell.type());
        const vector dX = xGp - xC;
        const scalar detJ =
            Foam::geometricJacobian::calcJacobianDetAtInteriorGaussPt
            (
                cell.type(),
                gaussEta[gpI],
                cellPoints
            );
        const scalar quadWeight = weights[gpI]*detJ;
        const Type qGp =
            evaluateLinearReconstruction
            (
                meanField[cellI],
                gradField[cellI],
                dX
            );

        forAll(rhs, dofI)
        {
            rhs[dofI] += qGp*(basis[gpI][dofI]*quadWeight);
        }
    }

    forAll(rhs, dofI)
    {
        if (mag(massDiag[dofI]) <= VSMALL)
        {
            FatalErrorInFunction
                << "Near-zero mass-matrix diagonal for field " << field.name()
                << " in cell " << cellI << ", dof " << dofI
                << ": " << massDiag[dofI] << exit(FatalError);
        }

        cellModes[dofI] = rhs[dofI]/massDiag[dofI];
    }

    field.dof().updateCellDof(cellI);
    field.gaussFields()[cellI].interpolateFromDof();
}


void convertConservativeToDG
(
    dgField<scalar>& rho,
    dgField<vector>& rhoU,
    dgField<scalar>& E,
    const volScalarField& pFoam,
    const volScalarField& TFoam,
    const volVectorField& UFoam,
    const dgGeomMesh& dgMesh,
    const eqnOfState& eos,
    const energy& energyModel,
    const bool heIsInternalEnergy
)
{
    volScalarField rhoFoam
    (
        IOobject
        (
            "rhoFoam",
            pFoam.time().timeName(),
            pFoam.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pFoam
    );

    volVectorField rhoUFoam
    (
        IOobject
        (
            "rhoUFoam",
            pFoam.time().timeName(),
            pFoam.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        rhoFoam*UFoam
    );

    volScalarField EFoam
    (
        IOobject
        (
            "EFoam",
            pFoam.time().timeName(),
            pFoam.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pFoam
    );

    forAll(rhoFoam, cellI)
    {
        const scalar rhoValue = eos.calcRhoFromPT(pFoam[cellI], TFoam[cellI]);
        const scalar TValue = TFoam[cellI];
        const vector UValue = UFoam[cellI];
        const scalar heValue =
            calcStoredThermoEnergy
            (
                rhoValue,
                TValue,
                eos,
                energyModel,
                heIsInternalEnergy
            );

        rhoFoam[cellI] = rhoValue;
        rhoUFoam[cellI] = rhoValue*UValue;
        EFoam[cellI] = rhoValue*(heValue + 0.5*magSqr(UValue));
    }

    buildBoundaryConservativeFields
    (
        rhoFoam,
        rhoUFoam,
        EFoam,
        pFoam,
        TFoam,
        UFoam,
        eos,
        energyModel,
        heIsInternalEnergy
    );

    tmp<volVectorField> tGradRho = fvc::grad(rhoFoam);
    tmp<volTensorField> tGradRhoU = fvc::grad(rhoUFoam);
    tmp<volVectorField> tGradE = fvc::grad(EFoam);

    const volVectorField& gradRho = tGradRho();
    const volTensorField& gradRhoU = tGradRhoU();
    const volVectorField& gradE = tGradE();

    forAll(rhoFoam, cellI)
    {
        projectLinearCellField(rho, cellI, rhoFoam, gradRho, dgMesh);
        projectLinearCellField(rhoU, cellI, rhoUFoam, gradRhoU, dgMesh);
        projectLinearCellField(E, cellI, EFoam, gradE, dgMesh);
    }

}

}
}


using namespace Foam;

int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, true);

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    if (timeDirs.size() != 1)
    {
        FatalErrorInFunction
            << "Please select exactly one target time with -time <time> "
            << "or -latestTime. Selected " << timeDirs.size()
            << " times." << exit(FatalError);
    }

    dgGeomMesh dgMesh(mesh);

    Info<< "Create DG geometric mesh with polynomial order "
        << dgMesh.pOrder() << nl << endl;

    runTime.setTime(timeDirs[0], 0);

    Info<< "Converting finite-volume fields at time "
        << runTime.timeName() << nl << endl;

    convertInputFoamPatchTypes(runTime);

    volScalarField pFoam
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );

    volScalarField TFoam
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );

    volVectorField UFoam
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );

    dgField<scalar> rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dgMesh,
        true
    );

    dgField<vector> rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dgMesh,
        true
    );

    dgField<scalar> E
    (
        IOobject
        (
            "E",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dgMesh,
        true
    );

    IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& dgThermoDict = thermoDict.subDict("dgThermo");
    const dictionary& mixDict = thermoDict.subDict("mixture");
    const word eosType = dgThermoDict.get<word>("equationOfState");
    const word thermoType = dgThermoDict.get<word>("thermo");
    const word energyType = dgThermoDict.get<word>("energy");

    autoPtr<eqnOfState> eos = eqnOfState::New(eosType, mixDict, dgMesh);
    autoPtr<thermoLaw> thermo = thermoLaw::New(thermoType, mixDict, dgMesh, eos());
    autoPtr<energy> energyModel = energy::New(energyType, mixDict, dgMesh, thermo());

    if (!(eos().isIdealGas() || eos().isThermalPerfectGas() || eos().isRealGas()))
    {
        FatalErrorInFunction
            << "convertFoamResultsToDG requires an ideal-gas, thermal-perfect, "
            << "or supported real-gas EOS, but got " << eos().type()
            << exit(FatalError);
    }

    if (!energyModel().heIsInternalEnergy() && !energyModel().heIsEnthalpy())
    {
        FatalErrorInFunction
            << "Unsupported energy model type " << energyModel().type()
            << exit(FatalError);
    }

    volScalarField rhoSeed
    (
        IOobject
        (
            "rhoSeed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pFoam
    );

    volVectorField rhoUSeed
    (
        IOobject
        (
            "rhoUSeed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        UFoam
    );

    volScalarField ESeed
    (
        IOobject
        (
            "ESeed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pFoam
    );

    forAll(pFoam, cellI)
    {
        const scalar rhoValue = eos().calcRhoFromPT(pFoam[cellI], TFoam[cellI]);
        const scalar heValue =
            calcStoredThermoEnergy
            (
                rhoValue,
                TFoam[cellI],
                eos(),
                energyModel(),
                energyModel().heIsInternalEnergy()
            );

        rhoSeed[cellI] = rhoValue;
        rhoUSeed[cellI] = rhoValue*UFoam[cellI];
        ESeed[cellI] = rhoValue*(heValue + 0.5*magSqr(UFoam[cellI]));
    }

    initialiseConservativeSeed(rho, rhoU, E, rhoSeed, ESeed, rhoUSeed);

    convertConservativeToDG
    (
        rho,
        rhoU,
        E,
        pFoam,
        TFoam,
        UFoam,
        dgMesh,
        eos(),
        energyModel(),
        energyModel().heIsInternalEnergy()
    );

    rho.write();
    rhoU.write();
    E.write();

    convertWrittenSymmetryPlanePatchTypes(runTime, dgMesh);

    Info<< "Wrote DG conservative restart fields at time "
        << runTime.timeName() << nl
        << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
