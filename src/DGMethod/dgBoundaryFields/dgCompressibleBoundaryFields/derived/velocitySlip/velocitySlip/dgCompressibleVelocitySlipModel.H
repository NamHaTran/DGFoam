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

#ifndef Foam_dgCompressibleVelocitySlipModel_H
#define Foam_dgCompressibleVelocitySlipModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "dgCompressibleBoundaryState.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "word.H"

namespace Foam
{

class dgGeomMesh;
class dgThermoConservative;
class fvPatch;

/*---------------------------------------------------------------------------*\
               Class dgCompressibleVelocitySlipModel Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleVelocitySlipModel
 * \brief Run-time selectable base class for wall velocity-slip models.
 *
 * \details
 * Derived models compute the effective boundary velocity at a boundary Gauss
 * point from the interior state and prescribed wall velocity.
 */
class dgCompressibleVelocitySlipModel
{
protected:

    /// Boundary patch using this slip model.
    const fvPatch& patch_;

    /// DG geometry mesh.
    const dgGeomMesh& dgMesh_;

    /// Conservative thermodynamic model and transport laws.
    const dgThermoConservative& thermo_;

    /// Model dictionary retained for diagnostics and derived classes.
    const dictionary dict_;

public:

    /// Run-time type name for the abstract model family.
    TypeName("velocitySlip");

    declareRunTimeSelectionTable
    (
        autoPtr,
        dgCompressibleVelocitySlipModel,
        dictionary,
        (
            const fvPatch& patch,
            const dgGeomMesh& dgMesh,
            const dgThermoConservative& thermo,
            const dictionary& dict
        ),
        (patch, dgMesh, thermo, dict)
    );

    /// Construct from a boundary patch, DG mesh, thermo model and dictionary.
    dgCompressibleVelocitySlipModel
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /// Destructor.
    virtual ~dgCompressibleVelocitySlipModel() = default;

    /// Create the named velocity-slip model from the dictionary table.
    static autoPtr<dgCompressibleVelocitySlipModel> New
    (
        const word& modelName,
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /**
     * \brief Return the effective boundary velocity.
     *
     * \param state Interior/exterior state at one boundary Gauss point.
     * \param Uwall Prescribed wall velocity.
     * \return Velocity used by the boundary condition.
     */
    virtual vector USlip
    (
        const dgCompressibleBoundaryState& state,
        const vector& Uwall
    ) const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
