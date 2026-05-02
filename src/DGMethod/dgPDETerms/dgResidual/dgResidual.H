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
    Foam::dgResidual

Description
    Cell-local residual assembler for a list of DG PDE terms.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgResidual_H
#define Foam_dgResidual_H

#include "FixedList.H"
#include "List.H"
#include "UPtrList.H"
#include "dgBasisField.H"
#include "dgGeomMesh.H"
#include "dgIntegralMath.H"
#include "dgPDETerm.H"
#include "error.H"
#include "pTraits.H"

#include <type_traits>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dgResidual Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgResidual
 * \brief Assemble one cell residual from multiple flux/source PDE terms.
 *
 * \tparam UType Type of the conservative variable and residual entry.
 * \tparam FType Type of the divergence flux associated with \c UType.
 *
 * \details
 * This class centralizes the residual integration that used to be repeated in
 * each individual PDE term.  For one equation on one cell, all convective,
 * diffusive, artificial-diffusive, and source terms are registered with
 * \c addTerm().  \c calcResidual() then performs one pass over the degrees of
 * freedom and Gauss points, summing all compatible term contributions before
 * applying the quadrature weights.
 *
 * For flux terms, the sign convention is
 *
 * \f[
 *     R_i^{F}
 *       =
 *       \int_K \nabla\phi_i \cdot F\,dV
 *       -
 *       \int_{\partial K} \phi_i F^\ast_n\,dS .
 * \f]
 *
 * For source terms, the contribution is
 *
 * \f[
 *     R_i^{S}
 *       =
 *       \int_K \phi_i S\,dV .
 * \f]
 *
 * Thus the final cell residual stored in \c R() is
 *
 * \f[
 *     R_i = R_i^F + R_i^S .
 * \f]
 *
 * Source terms are deliberately not added to the flux sum: they are multiplied
 * by \f$\phi_i\f$, while divergence fluxes are multiplied by
 * \f$\nabla\phi_i\f$ in the volume part and by \f$\phi_i\f$ only on faces.
 *
 * \note The assembler stores non-owning pointers to the terms.  Every term
 *       added through \c addTerm() must outlive the call to
 *       \c calcResidual().
 *
 * \note Storage uses fixed-size pointer lists to avoid dynamic allocation in
 *       the cell loop.  At most five PDE terms can be registered for one
 *       residual.
 */
template<class UType, class FType>
class dgResidual
{
    /** Maximum number of PDE terms supported by this compact assembler. */
    static constexpr label maxTerms_ = 5;

    /** Basis functions and derivatives for the owner cell. */
    const dgBasisField& basis_;

    /** Non-owning list of terms contributing to this cell residual. */
    FixedList<const dgPDETerm<UType, FType>*, maxTerms_> terms_;

    /** Number of active entries in \c terms_. */
    label nTerms_;

    /** Non-owning list of cell flux fields extracted from flux terms. */
    FixedList<const cellGaussField<FType>*, maxTerms_> fluxCellFields_;

    /** Non-owning list of face flux fields extracted from flux terms. */
    FixedList<const faceGaussField<FType>*, maxTerms_> fluxFaceFields_;

    /** Number of active entries in the flux-field lists. */
    label nFluxFields_;

    /** Non-owning list of source fields extracted from source terms. */
    FixedList<const cellGaussField<UType>*, maxTerms_> sourceFields_;

    /** Number of active entries in \c sourceFields_. */
    label nSourceFields_;

    /** Residual values indexed by local degree of freedom. */
    List<UType> R_;

public:

    /**
     * \brief Construct an empty residual assembler for one cell.
     *
     * Terms can be appended later with \c addTerm().
     */
    dgResidual(const dgGeomMesh& dgMesh, const label cellID)
    :
        basis_(dgMesh.basisField(cellID)),
        terms_(),
        nTerms_(0),
        fluxCellFields_(),
        fluxFaceFields_(),
        nFluxFields_(0),
        sourceFields_(),
        nSourceFields_(0),
        R_(basis_.nDof(), pTraits<UType>::zero)
    {}

    /**
     * \brief Construct a residual assembler from a non-owning pointer list.
     */
    dgResidual
    (
        const dgGeomMesh& dgMesh,
        const label cellID,
        const UPtrList<dgPDETerm<UType, FType>>& terms
    )
    :
        dgResidual(dgMesh, cellID)
    {
        forAll(terms, termI)
        {
            addTerm(terms[termI]);
        }
    }

    /**
     * \brief Register one PDE term in this residual.
     *
     * The term is not copied or owned.  It must remain valid until after
     * \c calcResidual() has been called.
     */
    void addTerm(const dgPDETerm<UType, FType>& term)
    {
        if (nTerms_ >= maxTerms_)
        {
            FatalErrorInFunction
                << "dgResidual supports at most " << maxTerms_
                << " PDE terms per cell residual" << abort(FatalError);
        }

        terms_[nTerms_++] = &term;

        if (const GaussField<FType>* F = term.fluxFieldPtr())
        {
            fluxCellFields_[nFluxFields_] = &F->cellField();
            fluxFaceFields_[nFluxFields_] = &F->faceField();
            ++nFluxFields_;
        }

        if (const cellGaussField<UType>* S = term.sourceFieldPtr())
        {
            sourceFields_[nSourceFields_++] = S;
        }
    }

    /** Return assembled residual values. */
    const List<UType>& R() const { return R_; }

    /** Return mutable assembled residual values. */
    List<UType>& R() { return R_; }

    /**
     * \brief Assemble the residual from all registered terms.
     *
     * The method clears \c R_ and then, for every basis function:
     * - sums all fluxes at each cell Gauss point before adding the volume
     *   contribution;
     * - sums all face numerical fluxes at each face Gauss point before adding
     *   the surface contribution;
     * - sums all source terms at each cell Gauss point before adding the source
     *   contribution.
     *
     * The final update keeps the DG sign convention explicit:
     *
     * \code
     * R_[dofI] += volTerm;
     * R_[dofI] -= surfTerm;
     * R_[dofI] += sourceTerm;
     * \endcode
     */
    void calcResidual()
    {
        forAll(R_, dofI)
        {
            R_[dofI] = pTraits<UType>::zero;
        }

        for (label dofI = 0; dofI < basis_.nDof(); ++dofI)
        {
            UType volTerm = pTraits<UType>::zero;
            UType surfTerm = pTraits<UType>::zero;
            UType sourceTerm = pTraits<UType>::zero;

            const GaussField<scalar>& phi = basis_.getBasis(dofI);

            if (dofI != 0)
            {
                const cellGaussField<vector>& gradPhi =
                    basis_.getDBasis(dofI).cellField();

                if (nFluxFields_ > 0)
                {
                    if constexpr (std::is_same_v<FType, tensor>)
                    {
                        volTerm = Foam::volumeIntegral
                        (
                            gradPhi,
                            fluxCellFields_,
                            nFluxFields_,
                            static_cast<const tensor*>(nullptr)
                        );
                    }
                    else
                    {
                        volTerm = Foam::volumeIntegral
                        (
                            gradPhi,
                            fluxCellFields_,
                            nFluxFields_
                        );
                    }
                }
            }

            if (nFluxFields_ > 0)
            {
                if constexpr (std::is_same_v<FType, tensor>)
                {
                    surfTerm = Foam::surfaceIntegral
                    (
                        phi.faceField(),
                        fluxFaceFields_,
                        nFluxFields_,
                        static_cast<const tensor*>(nullptr)
                    );
                }
                else
                {
                    surfTerm = Foam::surfaceIntegral
                    (
                        phi.faceField(),
                        fluxFaceFields_,
                        nFluxFields_
                    );
                }
            }

            if (nSourceFields_ > 0)
            {
                sourceTerm =
                    Foam::sourceIntegral<UType>
                    (
                        phi.cellField(),
                        sourceFields_,
                        nSourceFields_
                    );
            }

            R_[dofI] += volTerm;
            R_[dofI] -= surfTerm;
            R_[dofI] += sourceTerm;
        }
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
