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
    Foam::dgGeneralPDETerm

Description
    Generic DG PDE term used for residual assembly.

SourceFiles
    dgGeneralPDETerm.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeneralPDETerm_H
#define Foam_dgGeneralPDETerm_H

#include "word.H"
#include "List.H"
#include "label.H"
#include "pTraits.H"

#include "dgGeomMesh.H"
#include "dgBasisField.H"
#include "GaussField.H"
#include "dgFluxSolverManager.H"
#include "dgMath.H"

#include <type_traits>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dgGeneralPDETerm Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgGeneralPDETerm
 * \brief Generic DG PDE term for residual assembly
 *
 * \tparam UType Type of solution variable (scalar or vector)
 * \tparam FType Type of flux (vector or tensor)
 *
 * \details
 * This class represents a single PDE operator contribution to the DG
 * residual of one cell.
 *
 * The residual is assembled in the form:
 * \f[
 * R_i
 * =
 * - \int_{\Omega}
 *     \nabla \phi_i : \boldsymbol{F}\,\mathrm{d}V
 * + \int_{\partial\Omega}
 *     \phi_i (\boldsymbol{F}^{\ast}\cdot\boldsymbol{n})\,\mathrm{d}S
 * \f]
 *
 * Supported type combinations:
 * - UType = scalar, FType = vector
 * - UType = vector, FType = tensor
 *
 * The numerical surface flux \f$\boldsymbol{F}^{\ast}\f$ is computed
 * through dgFluxSolverManager based on the PDE term name.
 *
 * This class:
 * - Operates on a single cell
 * - Owns no global data
 * - Is used repeatedly during residual assembly
 */
template<class UType, class FType>
class dgGeneralPDETerm
{
private:

    /** Name of PDE term used to select flux solver */
    word termName_;

    /** Equation type for flux interpretation */
    dgFluxSolver::equationType eqnType_;

    /** Owner cell index */
    label cellID_;

    /** Reference to DG geometric mesh */
    dgGeomMesh& dgMesh_;

    /** Basis functions associated with this cell */
    dgBasisField basis_;

    /** Solution field evaluated at Gauss points */
    const GaussField<UType>& U_;

    /** Flux field evaluated at Gauss points */
    GaussField<FType>& F_;

    /** Flux solver manager */
    dgFluxSolverManager& fluxSolver_;

    /** Number of degrees of freedom in this cell */
    label nDof_;

    /** Residual values per degree of freedom */
    List<UType> R_;

    /**
     * Compile-time validation of UType/FType combination
     *
     * Ensures only supported scalar/vector formulations are instantiated.
     */
    inline void checkTypes() const;

public:

    // Constructors --------------------------------------------------------

    /**
     * \brief Construct PDE term for a given cell
     *
     * \param termName   Name of PDE term (used for flux solver lookup)
     * \param eqnType    Equation type (mass, momentum, energy, scalar)
     * \param cellID     Owner cell index
     * \param mesh       DG geometric mesh
     * \param U          Solution field at Gauss points
     * \param F          Flux field at Gauss points
     * \param fluxSolver Flux solver manager
     */
    dgGeneralPDETerm
    (
        const word& termName,
        dgFluxSolver::equationType eqnType,
        const label cellID,
        dgGeomMesh& mesh,
        const GaussField<UType>& U,
        GaussField<FType>& F,
        dgFluxSolverManager& fluxSolver
    );

    //- Destructor
    ~dgGeneralPDETerm() = default;

    // Access --------------------------------------------------------------

    /** Return number of degrees of freedom */
    inline label nDof() const { return nDof_; }

    /** Return residual values (const access) */
    inline const List<UType>& R() const { return R_; }

    /** Return residual values (mutable access) */
    inline List<UType>& R() { return R_; }

    // Core operations -----------------------------------------------------

    /**
     * \brief Compute numerical surface flux
     *
     * Selects the appropriate dgFluxSolver via dgFluxSolverManager
     * and updates the face Gauss flux field.
     */
    void updateFlux();

    /**
     * \brief Assemble cell residual
     *
     * Computes both volume and surface contributions to the
     * residual vector for this PDE term.
     */
    void calcResidual();
};

} // End namespace Foam

#ifdef NoRepository
    #include "dgGeneralPDETermI.H"
#endif

#endif

// ************************************************************************* //
