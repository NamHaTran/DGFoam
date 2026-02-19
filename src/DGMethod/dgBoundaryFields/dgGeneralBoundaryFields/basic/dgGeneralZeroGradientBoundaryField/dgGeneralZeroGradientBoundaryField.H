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

Class
    Foam::dgGeneralZeroGradientBoundaryField

Description
    Zero-gradient boundary condition for DG fields.

SourceFiles
    dgGeneralZeroGradientBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralZeroGradientBoundaryField_H
#define dgGeneralZeroGradientBoundaryField_H

#include "dgGeneralBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class dgGeneralZeroGradientBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeneralZeroGradientBoundaryField
 * \brief Zero-gradient (Neumann-type) boundary condition for DG fields
 *
 * \details
 * Enforces a zero normal gradient at face Gauss points by copying the
 * interior (minus-side) state to the ghost (plus-side) state.
 *
 * Ghost-state relation:
 * \f[
 *   U^{+} = U^{-}
 * \f]
 *
 * The gradient is assumed to be continuous across the boundary:
 * \f[
 *   \nabla U^{+} = \nabla U^{-}
 * \f]
 *
 * No dictionary parameters are required; only interior values are used.
 */
template<class Type>
class dgGeneralZeroGradientBoundaryField
:
    public dgGeneralBoundaryField<Type>
{
public:

    //- Runtime type information
    TypeName("generalZeroGradient");

    /**
     * \brief Construct zero-gradient boundary field
     *
     * \param patch  Finite-volume patch
     * \param dgMesh DG geometric mesh
     * \param dict   Boundary-condition dictionary (unused)
     */
    dgGeneralZeroGradientBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgGeneralZeroGradientBoundaryField() = default;
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgGeneralZeroGradientBoundaryFieldI.H"

#endif

// ************************************************************************* //
