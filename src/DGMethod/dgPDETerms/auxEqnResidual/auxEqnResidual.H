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

Class
    Foam::auxEqnResidual

Description
    Cell-local residual assembler for LDG auxiliary gradient-projection terms.

\*---------------------------------------------------------------------------*/

#ifndef Foam_auxEqnResidual_H
#define Foam_auxEqnResidual_H

#include "FixedList.H"
#include "List.H"
#include "UPtrList.H"
#include "dgBasisField.H"
#include "dgGradProjectionTerm.H"
#include "dgIntegralMath.H"
#include "error.H"
#include "pTraits.H"
#include "products.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class auxEqnResidual Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::auxEqnResidual
 * \brief Assemble auxiliary projection residuals for several state fields.
 *
 * \tparam UType Type of the projected conservative state \f$C\f$.
 *               Supported by \c dgGradProjectionTerm: \c scalar and \c vector.
 *
 * \details
 * This class mirrors the role of \c dgResidual for the auxiliary LDG block,
 * but each registered \c dgGradProjectionTerm produces a separate residual.
 * For every term \f$C_m\f$, it assembles
 *
 * \f[
 *     R_{m,i}
 *       =
 *      -\int_K C_m \nabla\phi_i\,dV
 *      +\int_{\partial K}\phi_i C_m^\ast n\,dS .
 * \f]
 *
 * Terms with the same \c UType share the same basis and quadrature loops.
 * For example, \c rho and \c E are both scalar projection terms and can be
 * assembled together, while \c rhoU is a vector projection term and is
 * assembled by a separate \c auxEqnResidual<vector>.
 *
 * \note The class stores non-owning pointers to terms and fields.  Registered
 *       terms must outlive \c calcResidual().
 *
 * \note Storage is fixed at five terms to avoid dynamic allocation in the
 *       cell loop.
 */
template<class UType>
class auxEqnResidual
{
public:

    /** Type of the projected state. */
    using state_type = UType;

    /** Type of the auxiliary residual entry. */
    using residual_type = typename outerProduct<UType, vector>::type;

private:

    /** Maximum number of auxiliary terms assembled together. */
    static constexpr label maxTerms_ = 5;

    /** Basis functions and derivatives for the owner cell. */
    const dgBasisField& basis_;

    /** Non-owning list of registered projection terms. */
    FixedList<const dgGradProjectionTerm<UType>*, maxTerms_> terms_;

    /** Non-owning list of cell state fields extracted from terms. */
    FixedList<const cellGaussField<UType>*, maxTerms_> cellFields_;

    /** Non-owning list of face state fields extracted from terms. */
    FixedList<const faceGaussField<UType>*, maxTerms_> faceFields_;

    /** Residual storage, indexed by [term][dof]. */
    FixedList<List<residual_type>, maxTerms_> R_;

    /** Number of active terms. */
    label nTerms_;

public:

    /** Construct an empty auxiliary residual assembler for one cell. */
    auxEqnResidual(const dgBasisField& basis)
    :
        basis_(basis),
        terms_(),
        cellFields_(),
        faceFields_(),
        R_(),
        nTerms_(0)
    {}

    /** Construct an assembler from a non-owning pointer list. */
    auxEqnResidual
    (
        const dgBasisField& basis,
        const UPtrList<dgGradProjectionTerm<UType>>& terms
    )
    :
        auxEqnResidual(basis)
    {
        forAll(terms, termI)
        {
            addTerm(terms[termI]);
        }
    }

    /**
     * \brief Register one auxiliary projection term.
     *
     * The returned residual index is the order in which the term was added.
     */
    label addTerm(const dgGradProjectionTerm<UType>& term)
    {
        if (nTerms_ >= maxTerms_)
        {
            FatalErrorInFunction
                << "auxEqnResidual supports at most " << maxTerms_
                << " terms per cell residual" << abort(FatalError);
        }

        const label termI = nTerms_++;

        terms_[termI] = &term;
        cellFields_[termI] = &term.stateField().cellField();
        faceFields_[termI] = &term.stateField().faceField();
        R_[termI].setSize(basis_.nDof(), pTraits<residual_type>::zero);

        return termI;
    }

    /** Return the number of active auxiliary terms. */
    label nTerms() const { return nTerms_; }

    /** Return assembled residual values for one registered term. */
    const List<residual_type>& R(const label termI) const
    {
        return R_[termI];
    }

    /** Return mutable residual values for one registered term. */
    List<residual_type>& R(const label termI)
    {
        return R_[termI];
    }

    /**
     * \brief Assemble residuals for all registered terms.
     *
     * The quadrature work is delegated to \c dgIntegralMath helpers:
     * \c gradProjectionVolumeIntegrals() and
     * \c gradProjectionCentralSurfaceIntegrals().
     */
    void calcResidual()
    {
        for (label termI = 0; termI < nTerms_; ++termI)
        {
            forAll(R_[termI], dofI)
            {
                R_[termI][dofI] = pTraits<residual_type>::zero;
            }
        }

        FixedList<residual_type, maxTerms_> volTerms;
        FixedList<residual_type, maxTerms_> surfTerms;

        for (label dofI = 0; dofI < basis_.nDof(); ++dofI)
        {
            for (label termI = 0; termI < nTerms_; ++termI)
            {
                volTerms[termI] = pTraits<residual_type>::zero;
            }

            if (dofI != 0)
            {
                Foam::gradProjectionVolumeIntegrals<UType>
                (
                    basis_.getDBasis(dofI).cellField(),
                    cellFields_,
                    nTerms_,
                    volTerms
                );
            }

            Foam::gradProjectionCentralSurfaceIntegrals<UType>
            (
                basis_.getBasis(dofI).faceField(),
                faceFields_,
                nTerms_,
                surfTerms
            );

            for (label termI = 0; termI < nTerms_; ++termI)
            {
                R_[termI][dofI] -= volTerms[termI];
                R_[termI][dofI] += surfTerms[termI];
            }
        }
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
