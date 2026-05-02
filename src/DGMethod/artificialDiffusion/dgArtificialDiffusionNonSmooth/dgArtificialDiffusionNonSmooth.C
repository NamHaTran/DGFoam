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

#include "dgArtificialDiffusionNonSmooth.H"
#include "basisFunctions.H"
#include "error.H"
#include "mathematicalConstants.H"

#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(dgArtificialDiffusionNonSmooth, 0);

IOobject dgArtificialDiffusionNonSmooth::createIOobject
(
    const dgGeomMesh& mesh
) const
{
    IOobject io
    (
        "dgOptions",
        mesh.getFvMesh().time().constant(),
        mesh.getFvMesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(false))
    {
        io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
    }
    else
    {
        io.readOpt(IOobject::NO_READ);
    }

    return io;
}


dgArtificialDiffusionNonSmooth::dgArtificialDiffusionNonSmooth
(
    const dgGeomMesh& mesh,
    const dgField<scalar>& rho,
    const dgField<vector>& U,
    const dgField<scalar>& a
)
:
    IOdictionary(createIOobject(mesh)),
    mesh_(mesh),
    rho_(rho),
    U_(U),
    a_(a),
    enabled_(false),
    type_("NonSmooth"),
    muAVName_("muAV"),
    mu0_(1.0),
    Skappa_(-1.0),
    Kappa_(0.25),
    sensorOffset_(1),
    diffusiveCFL_(0.25),
    report_(false),
    muAV_
    (
        IOobject
        (
            muAVName_,
            mesh.getFvMesh().time().timeName(),
            mesh.getFvMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        true
    ),
    tetNDof_(Zero),
    hexNDof_(Zero),
    prismNDof_(Zero),
    pyramidNDof_(Zero),
    tetHighModeDof_(),
    hexHighModeDof_(),
    prismHighModeDof_(),
    pyramidHighModeDof_()
{
    readCoeffs();

    tetHighModeDof_ =
        highOrderModeIndices
        (
            mesh_.pOrder(),
            dgCellType::TET,
            sensorOffset_,
            tetNDof_
        );
    hexHighModeDof_ =
        highOrderModeIndices
        (
            mesh_.pOrder(),
            dgCellType::HEX,
            sensorOffset_,
            hexNDof_
        );
    prismHighModeDof_ =
        highOrderModeIndices
        (
            mesh_.pOrder(),
            dgCellType::PRISM,
            sensorOffset_,
            prismNDof_
        );
    pyramidHighModeDof_ =
        highOrderModeIndices
        (
            mesh_.pOrder(),
            dgCellType::PYRAMID,
            sensorOffset_,
            pyramidNDof_
        );

    update();
}


void dgArtificialDiffusionNonSmooth::readCoeffs()
{
    if (readOpt() == IOobject::NO_READ || !found("artificialDiffusion"))
    {
        enabled_ = false;
        return;
    }

    const dictionary& coeffs = subDict("artificialDiffusion");

    enabled_ = coeffs.lookupOrDefault<bool>("enabled", false);
    type_ = coeffs.lookupOrDefault<word>("type", word("NonSmooth"));

    if (enabled_ && type_ != "NonSmooth")
    {
        FatalIOErrorInFunction(coeffs)
            << "Unsupported artificialDiffusion type '" << type_
            << "'. Only type 'NonSmooth' is implemented."
            << nl << exit(FatalIOError);
    }

    mu0_ = coeffs.lookupOrDefault<scalar>("mu0", mu0_);
    Skappa_ = coeffs.lookupOrDefault<scalar>("Skappa", Skappa_);
    Kappa_ = coeffs.lookupOrDefault<scalar>("Kappa", Kappa_);
    sensorOffset_ = coeffs.lookupOrDefault<label>("sensorOffset", sensorOffset_);
    diffusiveCFL_ =
        coeffs.lookupOrDefault<scalar>("diffusiveCFL", diffusiveCFL_);
    report_ = coeffs.lookupOrDefault<bool>("report", report_);

    if (sensorOffset_ < 1)
    {
        FatalIOErrorInFunction(coeffs)
            << "artificialDiffusion/sensorOffset must be >= 1."
            << nl << exit(FatalIOError);
    }

    if (Kappa_ <= VSMALL)
    {
        FatalIOErrorInFunction(coeffs)
            << "artificialDiffusion/Kappa must be positive."
            << nl << exit(FatalIOError);
    }

    if (diffusiveCFL_ <= SMALL)
    {
        FatalIOErrorInFunction(coeffs)
            << "artificialDiffusion/diffusiveCFL must be positive."
            << nl << exit(FatalIOError);
    }
}


const List<label>& dgArtificialDiffusionNonSmooth::highModeDof
(
    const dgCellType cellType
) const
{
    switch (cellType)
    {
        case dgCellType::TET:
            return tetHighModeDof_;

        case dgCellType::HEX:
            return hexHighModeDof_;

        case dgCellType::PRISM:
            return prismHighModeDof_;

        case dgCellType::PYRAMID:
            return pyramidHighModeDof_;

        default:
            FatalErrorInFunction
                << "Unsupported cell type " << cellType
                << " when accessing high-order modal indices."
                << abort(FatalError);
    }

    return tetHighModeDof_;
}


label dgArtificialDiffusionNonSmooth::expectedNDof
(
    const dgCellType cellType
) const
{
    switch (cellType)
    {
        case dgCellType::TET:
            return tetNDof_;

        case dgCellType::HEX:
            return hexNDof_;

        case dgCellType::PRISM:
            return prismNDof_;

        case dgCellType::PYRAMID:
            return pyramidNDof_;

        default:
            FatalErrorInFunction
                << "Unsupported cell type " << cellType
                << " when accessing modal size."
                << abort(FatalError);
    }

    return tetNDof_;
}


scalar dgArtificialDiffusionNonSmooth::densitySensor(const label cellID) const
{
    const dgGeomCell& cell = *mesh_.cells()[cellID];
    const dgCellType cellType = cell.type();
    const List<label>& highModes = highModeDof(cellType);
    const label expectedDof = expectedNDof(cellType);

    const GaussField<scalar>& rhoGF = rho_.gaussFields()[cellID];
    const cellGaussField<scalar>& rhoCell = rhoGF.cellField();
    const dofField<scalar>* rhoDofPtr = rhoGF.dofFieldPtr();

    if (!rhoDofPtr)
    {
        FatalErrorInFunction
            << "dgArtificialDiffusionNonSmooth requires density field '"
            << rho_.name() << "' to own modal DoFs."
            << abort(FatalError);
    }

    const List<scalar>& rhoDof = (*rhoDofPtr)[cellID].dof();

    if (rhoDof.size() != expectedDof)
    {
        FatalErrorInFunction
            << "Modal indexing mismatch for cell " << cellID
            << ": expected " << expectedDof << " modes but field '"
            << rho_.name() << "' stores " << rhoDof.size() << '.'
            << abort(FatalError);
    }

    const List<List<scalar>>& basis = cell.basis();
    const List<scalar>& weights = cell.weights();
    const List<scalar>& J3D = cell.J3D();

    scalar numerator = Zero;
    scalar denominator = Zero;

    forAll(weights, gpI)
    {
        scalar highOrderVal = Zero;

        forAll(highModes, modeI)
        {
            const label dofI = highModes[modeI];
            highOrderVal += basis[gpI][dofI]*rhoDof[dofI];
        }

        const scalar w = weights[gpI]*J3D[gpI];
        numerator += w*sqr(highOrderVal);
        denominator += w*sqr(rhoCell[gpI]);
    }

    if (denominator <= VSMALL || numerator <= VSMALL)
    {
        return std::log10(VSMALL);
    }

    return std::log10
    (
        max
        (
            std::sqrt(max(numerator/denominator, scalar(0))),
            VSMALL
        )
    );
}


scalar dgArtificialDiffusionNonSmooth::activation
(
    const scalar logSensor,
    const label pOrder
) const
{
    const scalar order = max(scalar(pOrder), scalar(1));
    const scalar Skappa = Skappa_ - scalar(4.25)*std::log10(order);

    if (logSensor < Skappa - Kappa_)
    {
        return Zero;
    }

    if (logSensor > Skappa + Kappa_)
    {
        return scalar(1);
    }

    return scalar(0.5)
       *(
            scalar(1)
          + std::sin(constant::mathematical::pi*(logSensor - Skappa)/(2*Kappa_))
        );
}


scalar dgArtificialDiffusionNonSmooth::hOverP(const label cellID) const
{
    const dgGeomCell& cell = *mesh_.cells()[cellID];
    const scalar p = max(scalar(mesh_.pOrder()), scalar(1));

    return max(cell.cellSize(), VSMALL)/p;
}


scalar dgArtificialDiffusionNonSmooth::maxWaveSpeed(const label cellID) const
{
    const cellGaussField<vector>& UCell = U_.gaussFields()[cellID].cellField();
    const cellGaussField<scalar>& aCell = a_.gaussFields()[cellID].cellField();

    scalar lambda = Zero;

    forAll(UCell, gpI)
    {
        lambda = max(lambda, mag(UCell[gpI]) + aCell[gpI]);
    }

    return lambda;
}


void dgArtificialDiffusionNonSmooth::setCellMuAV
(
    const label cellID,
    const scalar mu
)
{
    List<scalar>& cellDof = muAV_.dof()[cellID].dof();

    if (cellDof.empty())
    {
        FatalErrorInFunction
            << "muAV field has no DoFs in cell " << cellID << '.'
            << abort(FatalError);
    }

    forAll(cellDof, dofI)
    {
        cellDof[dofI] = scalar(0);
    }

    cellDof[0] = mu;
    muAV_.dof().updateCellDof(cellID);
}


dgField<scalar>& dgArtificialDiffusionNonSmooth::updateMuAV()
{
    if (!enabled_)
    {
        muAV_ = scalar(0);
        return muAV_;
    }

    for (label cellI = 0; cellI < mesh_.nCells(); ++cellI)
    {
        const scalar logSensor = densitySensor(cellI);
        const scalar chi = activation(logSensor, mesh_.pOrder());
        const scalar mu = chi*mu0_*maxWaveSpeed(cellI)*hOverP(cellI);

        setCellMuAV(cellI, mu);

        /*
        if (report_ && chi > scalar(0))
        {
            Info<< "artificialDiffusion NonSmooth: cell " << cellI
                << " logSensor=" << logSensor
                << " activation=" << chi
                << " muAV=" << mu << nl;
        }
        */
    }

    muAV_.interpolateFromDof();

    muAV_.synch();

    return muAV_;
}


scalar dgArtificialDiffusionNonSmooth::diffusiveDeltaT(const label cellID) const
{
    if (!enabled_)
    {
        return GREAT;
    }

    const dgGeomCell& cell = *mesh_.cells()[cellID];
    const scalar h = cell.cellSize();

    if (h <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive characteristic cell size for cell " << cellID
            << abort(FatalError);
    }

    const cellGaussField<scalar>& muAVCell =
        muAV_.gaussFields()[cellID].cellField();

    scalar muAVMax = Zero;

    forAll(muAVCell, gpI)
    {
        muAVMax = max(muAVMax, mag(muAVCell[gpI]));
    }

    if (muAVMax <= SMALL)
    {
        return GREAT;
    }

    const scalar pScale = sqr(scalar(mesh_.pOrder() + 1));
    const scalar diffusivePScale = sqr(pScale);

    return diffusiveCFL_*sqr(h)/(diffusivePScale*muAVMax);
}


void dgArtificialDiffusionNonSmooth::writeInfo(Ostream& os) const
{
    os  << "dgArtificialDiffusionNonSmooth: "
        << (enabled_ ? "enabled" : "disabled") << nl;

    if (enabled_)
    {
        os  << "  type         : " << type_ << nl
            << "  field        : " << muAV_.name() << nl
            << "  mu0          : " << mu0_ << nl
            << "  Skappa       : " << Skappa_ << nl
            << "  Kappa        : " << Kappa_ << nl
            << "  sensorOffset : " << sensorOffset_ << nl
            << "  diffusiveCFL : " << diffusiveCFL_ << nl;
    }
}

} // End namespace Foam

// ************************************************************************* //
