/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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

    /// Role of this patch in the coupled pair: \c inlet or \c outlet.
    word position_;

    /// Target outlet Knudsen number.
    scalar Kn_;

    /// Prescribed inlet Mach number.
    scalar MachIn_;

    /// Prescribed inlet temperature.
    scalar Tin_;

    /// Outlet reference temperature used to compute pressure from Kn.
    scalar TOut_;

    /// Inlet-to-outlet pressure ratio.
    scalar pressureRatio_;

    /// Reference length for Knudsen-number pressure evaluation.
    scalar refL_;

    /// Stored outlet pressure.
    scalar pOut_;

    /// Stored inlet pressure used when no outlet coupling is available.
    scalar pIn_;

    /// Optional prescribed inlet flow direction.
    vector direction_;

    /// Use the inward patch normal instead of a dictionary direction.
    bool usePatchNormalDirection_;

    /// Coupled inlet/outlet patch, resolved by the boundary manager.
    const dgCompressibleMachInletKnudsenOutletCouplingBoundaryField*
        coupledPatchPtr_;

    /// Compute outlet pressure from \c pOutValue or the target Knudsen number.
    scalar calcPout(const dictionary& coeffDict) const;

    /// Return the inlet velocity vector for face normal \p n.
    vector inletVelocity(const vector& n) const;

    /// Fill the conservative inlet ghost state.
    void updateInletState
    (
        const vector& n,
        scalar& rhoPlus,
        vector& rhoUPlus,
        scalar& EPlus
    ) const;

public:

    /// Run-time type name.
    TypeName("compressible::machInletKnudsenOutletCoupling");

    /// Construct from the patch, DG mesh, thermo model and dictionary.
    dgCompressibleMachInletKnudsenOutletCouplingBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /// Destructor.
    virtual ~dgCompressibleMachInletKnudsenOutletCouplingBoundaryField()
        = default;

    /// Return the configured patch role.
    inline const word& position() const { return position_; }

    /// Return the outlet pressure. Valid only for outlet-position patches.
    scalar pOut() const;

    /// Resolve and validate the inlet/outlet coupled patch pointer.
    virtual void resolveCoupling
    (
        const dgCompressibleBoundaryManager& manager
    ) override;

    /// Update the conservative ghost state for inlet or outlet operation.
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

    /// Update the primitive boundary value used for output/diagnostics.
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

    /// Ensure the boundary condition is attached to a generic patch.
    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
