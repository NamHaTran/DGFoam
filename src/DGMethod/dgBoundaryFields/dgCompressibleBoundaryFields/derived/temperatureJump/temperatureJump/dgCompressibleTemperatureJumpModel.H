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

#ifndef Foam_dgCompressibleTemperatureJumpModel_H
#define Foam_dgCompressibleTemperatureJumpModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "dgCompressibleBoundaryState.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "word.H"

namespace Foam
{

class dgGeomMesh;
class dgThermoConservative;
class fvPatch;

/*---------------------------------------------------------------------------*\
             Class dgCompressibleTemperatureJumpModel Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleTemperatureJumpModel
 * \brief Run-time selectable base class for wall temperature-jump models.
 *
 * \details
 * Derived models compute the effective boundary temperature at a boundary
 * Gauss point from the interior thermodynamic state and the prescribed wall
 * temperature.
 */
class dgCompressibleTemperatureJumpModel
{
protected:

    /// Boundary patch using this jump model.
    const fvPatch& patch_;

    /// DG geometry mesh.
    const dgGeomMesh& dgMesh_;

    /// Conservative thermodynamic model and transport laws.
    const dgThermoConservative& thermo_;

    /// Model dictionary retained for diagnostics and derived classes.
    const dictionary dict_;

public:

    /// Run-time type name for the abstract model family.
    TypeName("temperatureJump");

    declareRunTimeSelectionTable
    (
        autoPtr,
        dgCompressibleTemperatureJumpModel,
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
    dgCompressibleTemperatureJumpModel
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /// Destructor.
    virtual ~dgCompressibleTemperatureJumpModel() = default;

    /// Create the named temperature-jump model from the dictionary table.
    static autoPtr<dgCompressibleTemperatureJumpModel> New
    (
        const word& modelName,
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /**
     * \brief Return the effective boundary temperature.
     *
     * \param state Interior/exterior state at one boundary Gauss point.
     * \param Twall Prescribed wall temperature.
     * \return Temperature used by the boundary condition.
     */
    virtual scalar TJump
    (
        const dgCompressibleBoundaryState& state,
        const scalar Twall
    ) const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
