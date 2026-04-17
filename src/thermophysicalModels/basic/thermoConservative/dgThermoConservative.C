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

#include "dgThermoConservative.H"
#include "addToRunTimeSelectionTable.H"
#include "dgExpr.H"

namespace Foam
{
// * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dgThermoConservative, 0);
defineRunTimeSelectionTable(dgThermoConservative, dictionary);

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * //

Foam::dgThermoConservative::dgThermoConservative
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            "dgThermoConservative",                // object name
            mesh.getFvMesh().time().constant(),    // instance (constant/)
            mesh.getFvMesh().thisDb(),             // registry (like mesh)
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    dict_(dict),
    mesh_(mesh),

    // Null initialization. Derived class will allocate real models.
    eqnState_(nullptr),
    thermo_(nullptr),
    transport_(nullptr),
    energy_(nullptr),
    R_(Zero),
    he_
    (
        IOobject
        (
            "he",
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    a_
    (
        IOobject
        (
            "a", 
            mesh.getFvMesh().time().timeName(), 
            mesh.getFvMesh(), 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        false
    ),
    rho_
    (
        mesh_.getFvMesh().lookupObject<dgField<scalar>>("rho")
    ),
    rhoU_
    (
        mesh_.getFvMesh().lookupObject<dgField<vector>>("rhoU")
    ),
    E_
    (
        mesh_.getFvMesh().lookupObject<dgField<scalar>>("E")
    ),
    SRho_
    (
        mesh_.getFvMesh().lookupObject<dgField<vector>>("SRho")
    ),
    SRhoU_
    (
        mesh_.getFvMesh().lookupObject<dgField<tensor>>("SRhoU")
    ),
    SE_
    (
        mesh_.getFvMesh().lookupObject<dgField<vector>>("SE")
    ),
    U_
    (
        mesh_.getFvMesh().lookupObjectRef<dgField<vector>>("U")
    ),
    p_
    (
        mesh_.getFvMesh().lookupObjectRef<dgField<scalar>>("p")
    ),
    T_
    (
        mesh_.getFvMesh().lookupObjectRef<dgField<scalar>>("T")
    ),
    gradU_
    (
        mesh_.getFvMesh().lookupObject<dgField<tensor>>("gradU")
    ),
    gradP_
    (
        mesh_.getFvMesh().lookupObjectRef<dgField<vector>>("gradP")
    ),
    gradT_
    (
        mesh_.getFvMesh().lookupObjectRef<dgField<vector>>("gradT")
    ),
    TMean_(mesh_.nCells(), scalar(SMALL))
{
    // No model construction here.
    // Derived classes will call initModels() after base constructor completes.
}


// * * * * * * * * * * *  Factory Selector  * * * * * * * * * * * * //

Foam::autoPtr<Foam::dgThermoConservative> Foam::dgThermoConservative::New
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->find(name);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgThermoConservative type: " << name << nl
            << "Valid dgThermoConservative types: " << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(name, dict, mesh);
}

bool Foam::dgThermoConservative::writeData(Ostream& os) const
{
    return true;
}


void Foam::dgThermoConservative::synch()
{
    he_.synch();
    a_.synch();
}


void Foam::dgThermoConservative::synchGradient()
{
    gradP_.synch();
    gradT_.synch();
}


scalar Foam::dgThermoConservative::calcRawTemperatureFromRhoHe
(
    const scalar rho,
    const scalar he
) const
{
    if (heIsInternalEnergy() && eos().canCalcTFromRhoE())
    {
        return eos().calcTFromRhoE(rho, he);
    }

    return energyModel().calcTfromHe(he);
}


scalar Foam::dgThermoConservative::clampTemperature
(
    const scalar T,
    const scalar TMean
) const
{
    if (T >= 0)
    {
        return T;
    }

    if (TMean >= 0)
    {
        return TMean;
    }

    return scalar(SMALL);
}


scalar Foam::dgThermoConservative::clampTemperature
(
    const label cellID,
    const scalar T
) const
{
    return clampTemperature(T, TMean_[cellID]);
}


void Foam::dgThermoConservative::setTMean
(
    const label cellID,
    const scalar TMean
)
{
    TMean_[cellID] = clampTemperature(TMean, TMean_[cellID]);
}


scalar Foam::dgThermoConservative::calcTemperatureFromRhoHe
(
    const label cellID,
    const scalar rho,
    const scalar he
) const
{
    return clampTemperature(cellID, calcRawTemperatureFromRhoHe(rho, he));
}


scalar Foam::dgThermoConservative::calcPressureFromRhoHe
(
    const label cellID,
    const scalar rho,
    const scalar he
) const
{
    if (heIsInternalEnergy() && eos().canCalcPFromRhoE())
    {
        return eos().calcPFromRhoE(rho, he);
    }

    return eos().calcPFromRhoT(rho, calcTemperatureFromRhoHe(cellID, rho, he));
}


scalar Foam::dgThermoConservative::calcSpeedOfSoundFromRhoHe
(
    const label cellID,
    const scalar rho,
    const scalar he
) const
{
    if (heIsInternalEnergy() && eos().canCalcAFromRhoE())
    {
        return eos().calcAFromRhoE(rho, he);
    }

    const scalar T = calcTemperatureFromRhoHe(cellID, rho, he);
    const scalar Cp = thermo().calcCp(T);
    const scalar Cv = thermo().calcCv(T);
    const scalar gamma = thermo().calcGamma(Cp, Cv);

    return thermo().calcSpeedOfSound(T, gamma);
}


scalar Foam::dgThermoConservative::calcHeFromRhoT
(
    const scalar rho,
    const scalar T
) const
{
    if (heIsInternalEnergy() && eos().canCalcEFromRhoT())
    {
        return eos().calcEFromRhoT(rho, T);
    }

    return energyModel().calcHe(T);
}


tmp<GaussField<scalar>> Foam::dgThermoConservative::calcH
(
    const label cellID
) const
{
    tmp<GaussField<scalar>> th = GaussField<scalar>::New(cellID, &mesh_);
    GaussField<scalar>& h = th.ref();

    if (heIsEnthalpy())
    {
        h = he_.gaussFields()[cellID];
    }
    else
    {
        const GaussField<scalar>& he = he_.gaussFields()[cellID];
        const GaussField<scalar>& p = p_.gaussFields()[cellID];
        const GaussField<scalar>& rho = rho_.gaussFields()[cellID];
        dg::assign(h, dg::expr(he) + dg::expr(p)/dg::expr(rho));
    }

    return th;
}


void Foam::dgThermoConservative::calcH
(
    const label cellID,
    const label localFaceID,
    boundaryGaussField<scalar>& hMinus,
    boundaryGaussField<scalar>& hPlus
) const
{
    const faceGaussField<scalar>& he = he_.gaussFields()[cellID].faceField();
    const label nGauss = he_.gaussFields()[cellID].faceField().nGaussPerFace();

    if (hMinus.size() != nGauss)
    {
        hMinus = boundaryGaussField<scalar>(nGauss);
    }

    if (hPlus.size() != nGauss)
    {
        hPlus = boundaryGaussField<scalar>(nGauss);
    }

    if (heIsEnthalpy())
    {
        for (label gpI = 0; gpI < nGauss; ++gpI)
        {
            hMinus[gpI] = he.minusValueOnFace(localFaceID, gpI);
            hPlus[gpI] = he.plusValueOnFace(localFaceID, gpI);
        }
    }
    else
    {
        const faceGaussField<scalar>& p = p_.gaussFields()[cellID].faceField();
        const faceGaussField<scalar>& rho = rho_.gaussFields()[cellID].faceField();

        for (label gpI = 0; gpI < nGauss; ++gpI)
        {
            hMinus[gpI] =
                he.minusValueOnFace(localFaceID, gpI)
              + p.minusValueOnFace(localFaceID, gpI)
               /rho.minusValueOnFace(localFaceID, gpI);

            hPlus[gpI] =
                he.plusValueOnFace(localFaceID, gpI)
              + p.plusValueOnFace(localFaceID, gpI)
               /rho.plusValueOnFace(localFaceID, gpI);
        }
    }
}


tmp<GaussField<scalar>> Foam::dgThermoConservative::calcE
(
    const label cellID
) const
{
    tmp<GaussField<scalar>> te = GaussField<scalar>::New(cellID, &mesh_);
    GaussField<scalar>& e = te.ref();

    if (heIsInternalEnergy())
    {
        e = he_.gaussFields()[cellID];
    }
    else
    {
        const GaussField<scalar>& he = he_.gaussFields()[cellID];
        const GaussField<scalar>& p = p_.gaussFields()[cellID];
        const GaussField<scalar>& rho = rho_.gaussFields()[cellID];
        dg::assign(e, dg::expr(he) - dg::expr(p)/dg::expr(rho));
    }

    return te;
}


void Foam::dgThermoConservative::calcE
(
    const label cellID,
    const label localFaceID,
    boundaryGaussField<scalar>& eMinus,
    boundaryGaussField<scalar>& ePlus
) const
{
    const faceGaussField<scalar>& he = he_.gaussFields()[cellID].faceField();
    const label nGauss = he_.gaussFields()[cellID].faceField().nGaussPerFace();

    if (eMinus.size() != nGauss)
    {
        eMinus = boundaryGaussField<scalar>(nGauss);
    }

    if (ePlus.size() != nGauss)
    {
        ePlus = boundaryGaussField<scalar>(nGauss);
    }

    if (heIsInternalEnergy())
    {
        for (label gpI = 0; gpI < nGauss; ++gpI)
        {
            eMinus[gpI] = he.minusValueOnFace(localFaceID, gpI);
            ePlus[gpI] = he.plusValueOnFace(localFaceID, gpI);
        }
    }
    else
    {
        const faceGaussField<scalar>& p = p_.gaussFields()[cellID].faceField();
        const faceGaussField<scalar>& rho = rho_.gaussFields()[cellID].faceField();

        for (label gpI = 0; gpI < nGauss; ++gpI)
        {
            eMinus[gpI] =
                he.minusValueOnFace(localFaceID, gpI)
              - p.minusValueOnFace(localFaceID, gpI)
               /rho.minusValueOnFace(localFaceID, gpI);

            ePlus[gpI] =
                he.plusValueOnFace(localFaceID, gpI)
              - p.plusValueOnFace(localFaceID, gpI)
               /rho.plusValueOnFace(localFaceID, gpI);
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
