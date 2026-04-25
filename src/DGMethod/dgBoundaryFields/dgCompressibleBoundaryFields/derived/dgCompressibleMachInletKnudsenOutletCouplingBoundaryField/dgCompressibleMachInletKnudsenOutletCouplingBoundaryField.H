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

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCompressibleMachInletKnudsenOutletCouplingBoundaryField_H
#define Foam_dgCompressibleMachInletKnudsenOutletCouplingBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

class dgCompressibleBoundaryManager;

/*---------------------------------------------------------------------------*\
 Class dgCompressibleMachInletKnudsenOutletCouplingBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleMachInletKnudsenOutletCouplingBoundaryField
 * \brief Coupled Mach inlet / Knudsen outlet boundary condition.
 *
 * \details
 * The outlet side stores a pressure computed from the target Knudsen number
 * and applies a zero-gradient conservative ghost state.  The inlet side can
 * couple to that outlet patch, read its stored outlet pressure, derive the
 * inlet pressure from a pressure ratio, and impose fixed p/T/U.
 */
class dgCompressibleMachInletKnudsenOutletCouplingBoundaryField
:
    public dgCompressibleBoundaryField
{
private:

    word position_;

    scalar Kn_;
    scalar MachIn_;
    scalar Tin_;
    scalar TOut_;
    scalar pressureRatio_;
    scalar refL_;

    scalar pOut_;
    scalar pIn_;

    vector direction_;
    bool usePatchNormalDirection_;

    const dgCompressibleMachInletKnudsenOutletCouplingBoundaryField*
        coupledPatchPtr_;

    scalar calcPout(const dictionary& coeffDict) const;

    vector inletVelocity(const vector& n) const;

    void updateInletState
    (
        const vector& n,
        scalar& rhoPlus,
        vector& rhoUPlus,
        scalar& EPlus
    ) const;

public:

    TypeName("compressible::machInletKnudsenOutletCoupling");

    dgCompressibleMachInletKnudsenOutletCouplingBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    virtual ~dgCompressibleMachInletKnudsenOutletCouplingBoundaryField()
        = default;

    inline const word& position() const { return position_; }

    scalar pOut() const;

    virtual void resolveCoupling
    (
        const dgCompressibleBoundaryManager& manager
    ) override;

    virtual void updateConservativeGhostState
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar rhoMinus,
        const vector& rhoUMinus,
        const scalar EMinus,
        scalar& rhoPlus,
        vector& rhoUPlus,
        scalar& EPlus
    ) const override;

    virtual void updatePrimitiveBCValue
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar rhoMinus,
        const vector& rhoUMinus,
        const scalar EMinus,
        scalar& rhoBC,
        vector& rhoUBC,
        scalar& EBC
    ) const override;

    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
