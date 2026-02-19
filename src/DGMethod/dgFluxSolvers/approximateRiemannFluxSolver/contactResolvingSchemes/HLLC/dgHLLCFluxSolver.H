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
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dgHLLCFluxSolver

Description
    HLLC approximate Riemann solver for DG methods.

SourceFiles
    dgHLLCFluxSolver.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_dgHLLCFluxSolver_H
#define Foam_dgHLLCFluxSolver_H

#include "dgFluxSolver.H"
#include "dgField.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{
    template<class Type> class faceGaussField;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dgHLLCFluxSolver Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgHLLCFluxSolver
 * \brief HLLC approximate Riemann solver for DG methods
 *
 * \details
 * Implements the HLLC (Harten–Lax–van Leer–Contact) numerical flux.
 *
 * HLLC extends HLLE by restoring the contact (and shear) wave through a
 * middle wave speed \f$S^*\f$ and two star states. This reduces numerical
 * dissipation and improves resolution of contact discontinuities.
 *
 * The local one-dimensional Riemann problem is solved along the face unit
 * normal \f$\boldsymbol{n}\f$.
 *
 * Flux selection by wave-speed ordering:
 * \f[
 * \boldsymbol{f} =
 * \begin{cases}
 * \boldsymbol{f}_L,        & S_L \ge 0, \\
 * \boldsymbol{f}_L^*,      & S_L \le 0 \le S^*, \\
 * \boldsymbol{f}_R^*,      & S^* \le 0 \le S_R, \\
 * \boldsymbol{f}_R,        & S_R \le 0
 * \end{cases}
 * \f]
 *
 * Cartesian flux handling:
 * - Scalar transport: returns \f$f\,\boldsymbol{n}\f$
 * - Momentum transport: returns vector flux directly
 */
class dgHLLCFluxSolver
:
    public dgFluxSolver
{
public:

    /** Wave-speed estimation methods */
    enum speedEstimateType
    {
        seDavis,        ///< Davis estimate (Un ± a)
        seRoeEinfeldt   ///< Roe–Einfeldt averaged estimate
    };

private:

    /** Pointer to density field (read-only, non-owning) */
    const dgField<scalar>& rho_;

    /** Pointer to velocity field (read-only, non-owning) */
    const dgField<vector>& U_;

    /** Pointer to pressure field (read-only, non-owning) */
    const dgField<scalar>& p_;

    /** Pointer to speed-of-sound field (read-only, non-owning) */
    const dgField<scalar>& a_;

    /** Pointer to enthalpy field (read-only, non-owning) */
    const dgField<scalar>& h_;

    /** Pointer to heat-capacity ratio field (read-only, non-owning) */
    const dgField<scalar>& gamma_;

    /** Selected wave-speed estimation method */
    speedEstimateType speedEst_;

    /** Cached left wave speeds per face and Gauss point */
    List<List<scalar>> SL_list_;

    /** Cached right wave speeds per face and Gauss point */
    List<List<scalar>> SR_list_;

    /** Cached contact wave speeds per face and Gauss point */
    List<List<scalar>> SStar_list_;

    /** Cached left star-state coefficients */
    List<List<scalar>> CL_list_;

    /** Cached right star-state coefficients */
    List<List<scalar>> CR_list_;

    /** Flags indicating whether intermediate states are computed */
    List<bool> isStateComputed_;

public:

    //- Runtime type information
    TypeName("HLLC");

    /**
     * \brief Construct from name, dictionary, and DG mesh
     *
     * Reads wave-speed estimation settings and allocates storage for
     * cached intermediate HLLC states.
     *
     * \param name Solver name
     * \param dict Flux-solver dictionary
     * \param mesh DG geometric mesh
     */
    dgHLLCFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //- Destructor
    virtual ~dgHLLCFluxSolver() = default;

    /**
     * \brief Compute intermediate wave speeds and star-state coefficients
     *
     * Computes \f$S_L\f$, \f$S_R\f$, \f$S^\*\f$ and star-state scaling
     * coefficients \f$C_L\f$, \f$C_R\f$ at a face Gauss point.
     *
     * \param cellID        DG cell index
     * \param localFaceID   Local face index within the cell
     * \param localGaussID  Local Gauss-point index on the face
     * \param n             Face unit normal (minus-to-plus)
     * \param SL            Left wave speed (output)
     * \param SR            Right wave speed (output)
     * \param SStar         Contact wave speed (output)
     * \param CL            Left star-state coefficient (output)
     * \param CR            Right star-state coefficient (output)
     */
    void calcIntermediateState
    (
        const label cellID,
        const label localFaceID,
        const label localGaussID,
        const vector& n,
        scalar& SL,
        scalar& SR,
        scalar& SStar,
        scalar& CL,
        scalar& CR
    );

    /**
     * \brief Compute numerical flux for scalar conserved variable
     *
     * \param cellID DG cell index
     * \param F      Flux field (output)
     * \param U      Scalar solution at face Gauss points
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;

    /**
     * \brief Compute numerical flux for vector conserved variable
     *
     * \param cellID DG cell index
     * \param F      Flux field (output)
     * \param U      Vector solution at face Gauss points
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;

    /**
     * \brief Reset cached intermediate HLLC states
     */
    virtual void reset() override;

    /**
     * \brief Read solver settings from dictionary
     *
     * \param dict Flux-solver dictionary
     */
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif
// ************************************************************************* //
