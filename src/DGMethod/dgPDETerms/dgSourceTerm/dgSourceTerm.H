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
    Foam::dgSourceTerm

Description
    Volume-source DG PDE term used to expose source fields.

SourceFiles
    dgSourceTermI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgSourceTerm_H
#define Foam_dgSourceTerm_H

#include "word.H"
#include "List.H"
#include "label.H"
#include "pTraits.H"

#include "dgBasisField.H"
#include "dgGeomMesh.H"
#include "GaussField.H"
#include "dgPDETerm.H"
#include "products.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dgSourceTerm Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgSourceTerm
 * \brief Cell-local DG volume source term.
 *
 * \tparam Type Type of source and target conservative variable.
 *
 * \details
 * This class stores a cell-local source field and exposes it to dgResidual.
 * It has no face contribution and no numerical flux. For a conservation law
 *
 * \f[
 *     \partial_t Q + \nabla\cdot F = S,
 * \f]
 *
 * add this source contribution directly to the matching residual.
 * If a model writes the source on the left-hand side, pass the negated source
 * field instead.
 */
template<class Type, class FType = typename outerProduct<Type, vector>::type>
class dgSourceTerm
:
    public dgPDETerm<Type, FType>
{
private:

    /** Name of PDE term, useful for diagnostics and symmetry with flux terms. */
    word termName_;

    /** DG mesh used to retrieve cell-local data. */
    const dgGeomMesh& dgMesh_;

    /** Active source field evaluated at cell-interior Gauss points. */
    const cellGaussField<Type>* S_;

public:

    using source_type = Type;

    dgSourceTerm
    (
        const word& termName,
        const dgGeomMesh& dgMesh
    );

    void add
    (
        const cellGaussField<Type>& S
    );

    void add
    (
        const GaussField<Type>& S
    );

    void clear();

    ~dgSourceTerm() = default;

    inline bool isSource() const override { return true; }

    inline const cellGaussField<Type>& sourceField() const { return *S_; }

    inline const cellGaussField<Type>* sourceFieldPtr() const override
    {
        return S_;
    }

    inline const word& termName() const { return termName_; }
};

} // End namespace Foam

#ifdef NoRepository
    #include "dgSourceTermI.H"
#endif

#endif

// ************************************************************************* //
