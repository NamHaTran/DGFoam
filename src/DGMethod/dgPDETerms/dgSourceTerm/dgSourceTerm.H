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
    Volume-source DG PDE term used for residual assembly.

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
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dgSourceTerm Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgSourceTerm
 * \brief Cell-local DG volume source term for residual assembly.
 *
 * \tparam Type Type of source and target conservative variable.
 *
 * \details
 * This class projects a source field \f$S\f$ onto the DG basis:
 *
 * \f[
 *     R_i = \int_{\Omega_K} \phi_i S\,\mathrm{d}V.
 * \f]
 *
 * It has no face contribution and no numerical flux. For a conservation law
 *
 * \f[
 *     \partial_t Q + \nabla\cdot F = S,
 * \f]
 *
 * add this residual directly to the matching time-discretization buffer.
 * If a model writes the source on the left-hand side, pass the negated source
 * field instead.
 */
template<class Type>
class dgSourceTerm
{
private:

    /** Name of PDE term, useful for diagnostics and symmetry with flux terms. */
    word termName_;

    /** Owner cell index. */
    label cellID_;

    /** Basis functions associated with this cell. */
    const dgBasisField& basis_;

    /** Source field evaluated at cell-interior Gauss points. */
    tmp<cellGaussField<Type>> S_;

    /** Number of degrees of freedom in this cell. */
    label nDof_;

    /** Residual values per degree of freedom. */
    List<Type> R_;

public:

    using source_type = Type;

    dgSourceTerm
    (
        const word& termName,
        const dgBasisField& basis,
        const cellGaussField<Type>& S
    );

    dgSourceTerm
    (
        const word& termName,
        const dgBasisField& basis,
        const tmp<cellGaussField<Type>>& tS
    );

    dgSourceTerm
    (
        const word& termName,
        const dgBasisField& basis,
        const GaussField<Type>& S
    );

    ~dgSourceTerm() = default;

    inline label nDof() const { return nDof_; }

    inline const List<Type>& R() const { return R_; }

    inline List<Type>& R() { return R_; }

    inline const cellGaussField<Type>& sourceField() const { return S_(); }

    inline const word& termName() const { return termName_; }

    /** \brief Assemble \f$R_i = \int \phi_i S\,dV\f$ on the owner cell. */
    void calcResidual();
};

} // End namespace Foam

#ifdef NoRepository
    #include "dgSourceTermI.H"
#endif

#endif

// ************************************************************************* //
