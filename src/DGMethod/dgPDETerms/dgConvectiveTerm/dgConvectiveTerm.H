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
    Foam::dgConvectiveTerm

Description
    Convective DG PDE term used to compute and expose fluxes.

SourceFiles
    dgConvectiveTerm.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgConvectiveTerm_H
#define Foam_dgConvectiveTerm_H

#include "word.H"
#include "List.H"
#include "label.H"
#include "pTraits.H"

#include "dgGeomMesh.H"
#include "dgBasisField.H"
#include "GaussField.H"
#include "dgFluxSolverManager.H"
#include "dgMath.H"
#include "dgPDETerm.H"

#include <type_traits>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dgConvectiveTerm Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgConvectiveTerm
 * \brief Convective DG PDE term for flux assembly.
 *
 * \tparam UType Type of solution variable (scalar or vector)
 * \tparam FType Type of flux (vector or tensor)
 *
 * \details
 * This class represents a single PDE operator contribution on one cell.
 * It updates the numerical face flux and exposes the Gauss flux field to
 * dgResidual, which performs the common residual integration.
 *
 * Supported type combinations:
 * - UType = scalar, FType = vector
 * - UType = vector, FType = tensor
 *
 * The numerical surface flux \f$\boldsymbol{F}^{\ast}\f$ is computed
 * through dgFluxSolverManager based on the PDE term name.
 *
 * This class owns no global data and is used repeatedly during residual
 * assembly.
 */
template<class UType, class FType>
class dgConvectiveTerm
:
    public dgPDETerm<UType, FType>
{
private:

    /** Name of PDE term used to select flux solver */
    word termName_;

    /** Equation type for flux interpretation */
    dgFluxSolver::equationType eqnType_;

    /** Flux type required by this PDE term */
    dgFluxSolver::fluxType fluxType_;

    /** Owner cell index */
    label cellID_;

    /** Basis functions associated with this cell */
    const dgBasisField& basis_;

    /** Solution field evaluated at Gauss points */
    const GaussField<UType>& U_;

    /** Flux field evaluated at Gauss points */
    tmp<GaussField<FType>> F_;

    /** Flux solver manager */
    dgFluxSolverManager& fluxSolver_;

    /** Flux solver selected for this term. */
    dgFluxSolver& termFluxSolver_;

    /**
     * Compile-time validation of UType/FType combination
     *
     * Ensures only supported scalar/vector formulations are instantiated.
     */
    inline void checkTypes() const;

public:

    using solution_type = UType;
    using flux_type = FType;

    // Constructors --------------------------------------------------------

    /**
     * \brief Construct PDE term for a given cell
     *
     * \param termName   Name of PDE term (used for flux solver lookup)
     * \param eqnType    Equation type (mass, momentum, energy, scalar)
     * \param basis      Basis field associated with the owner cell
     * \param U          Solution field at Gauss points
     * \param F          Flux field at Gauss points
     * \param fluxSolver Flux solver manager
     */
    dgConvectiveTerm
    (
        const word& termName,
        dgFluxSolver::equationType eqnType,
        dgFluxSolver::fluxType fluxType,
        const dgBasisField& basis,
        const GaussField<UType>& U,
        GaussField<FType>& F,
        dgFluxSolverManager& fluxSolver
    );

    /**
     * \brief Construct PDE term for a given cell, with Flux field is tmp
     *
     * \param termName   Name of PDE term (used for flux solver lookup)
     * \param eqnType    Equation type (mass, momentum, energy, scalar)
     * \param basis      Basis field associated with the owner cell
     * \param U          Solution field at Gauss points
     * \param F          Flux field at Gauss points
     * \param fluxSolver Flux solver manager
     */
    dgConvectiveTerm
    (
        const word& termName,
        dgFluxSolver::equationType eqnType,
        dgFluxSolver::fluxType fluxType,
        const dgBasisField& basis,
        const GaussField<UType>& U,
        const tmp<GaussField<FType>>& tF,
        dgFluxSolverManager& fluxSolver
    );

    //- Destructor
    ~dgConvectiveTerm() = default;

    // Access --------------------------------------------------------------

    inline bool isConvection() const override { return true; }

    /** Return Gauss flux field, including face numerical flux values */
    inline const GaussField<FType>& fluxField() const { return F_(); }

    inline const GaussField<FType>* fluxFieldPtr() const override
    {
        return &F_();
    }

    /** Return flux type required by this term */
    inline dgFluxSolver::fluxType fluxType() const { return fluxType_; }

    // Core operations -----------------------------------------------------

    /**
     * \brief Compute numerical surface flux
     *
     * Selects the appropriate dgFluxSolver via dgFluxSolverManager
     * and updates the face Gauss flux field.
     */
    void updateFlux();

};

} // End namespace Foam

#ifdef NoRepository
    #include "dgConvectiveTermI.H"
#endif

#endif

// ************************************************************************* //
