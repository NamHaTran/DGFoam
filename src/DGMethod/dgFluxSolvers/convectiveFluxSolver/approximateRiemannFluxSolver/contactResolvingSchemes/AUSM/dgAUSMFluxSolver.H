/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

Class
    Foam::dgAUSMFluxSolver

Description
    AUSM approximate Riemann solver for DG methods.

SourceFiles
    dgAUSMFluxSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgAUSMFluxSolver_H
#define Foam_dgAUSMFluxSolver_H

#include "dgFluxSolver.H"
#include "dgField.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{
    class dgThermoConservative;
    template<class Type> class faceGaussField;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dgAUSMFluxSolver Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgAUSMFluxSolver
 * \brief AUSM-family numerical flux for the compressible Euler equations
 *
 * \details
 * This solver evaluates the interface flux in a face-local orthonormal frame
 * whose first axis is aligned with the face normal. In that frame, the Euler
 * flux is split into:
 * - an advective part controlled by a split interface Mach number, and
 * - a pressure part reconstructed from split pressure polynomials.
 *
 * The implementation exposes a single runtime-selectable solver name,
 * \c AUSM, with a \c variant option to switch among supported coefficient
 * sets:
 * - \c AUSM0
 * - \c AUSM1
 * - \c AUSM2
 * - \c AUSM3
 *
 * Supported transport equations:
 * - massTransport
 * - momentumTransport
 * - energyTransport
 *
 * Generic scalarTransport is intentionally rejected because this solver is
 * formulated for the coupled Euler system rather than passive-scalar fluxing.
 *
 * Example dictionary usage:
 * \verbatim
 * fluxSchemes
 * {
 *     fluxSolvers
 *     {
 *         defaultConvective   AUSM;
 *         massConvTermFlux    AUSM;
 *         momentumConvTermFlux AUSM;
 *         energyConvTermFlux  AUSM;
 *     }
 *
 *     fluxSolversCoeffs
 *     {
 *         AUSMCoeffs
 *         {
 *             variant AUSM3;
 *             Mco     0.01;
 *         }
 *     }
 * }
 * \endverbatim
 */
class dgAUSMFluxSolver
:
    public dgFluxSolver
{
public:

    /** Available AUSM coefficient sets. */
    enum class variantType
    {
        ausm0,  //!< Baseline AUSM splitting without extra correction terms.
        ausm1,  //!< AUSM with smoother subsonic Mach/pressure splitting.
        ausm2,  //!< AUSM with pressure- and velocity-coupling corrections.
        ausm3   //!< AUSM2-style corrections with low-Mach scaling in Mp.
    };

protected:

    /** Thermodynamic package used to recover sound speed and primitive data. */
    const dgThermoConservative& thermo_;

    /** Conservative density field. */
    const dgField<scalar>& rho_;

    /** Conservative total energy field. */
    const dgField<scalar>& E_;

    /** Primitive velocity field. */
    const dgField<vector>& U_;

    /** Primitive pressure field. */
    const dgField<scalar>& p_;

    /** Selected AUSM variant. */
    variantType variant_;

    /** Polynomial coefficient used in the Mach-number split. */
    scalar beta_;

    /** Polynomial coefficient used in the pressure split. */
    scalar alpha_;

    /** Pressure-difference coefficient entering the Mach correction term. */
    scalar Kp_;

    /** Velocity-difference coefficient entering the pressure correction term. */
    scalar Ku_;

    /** Sensor coefficient used by the pressure-difference Mach correction. */
    scalar sigma_;

    /** Lower bound used by the AUSM3 low-Mach scaling. */
    scalar Mco_;

    /** Enable the pressure-difference contribution to the interface Mach number. */
    bool usePressureMachCorrection_;

    /** Enable the velocity-difference contribution to the interface pressure. */
    bool useVelocityPressureCorrection_;

    /**
     * \brief Positive/negative supersonic Mach split.
     * \param side Select plus/minus polynomial branch.
     * \param M    Local Mach number.
     * \return Supersonic split contribution.
     */
    scalar M1Function
    (
        const label side,
        const scalar M
    ) const;

    /**
     * \brief Positive/negative subsonic Mach polynomial.
     * \param side Select plus/minus polynomial branch.
     * \param M    Local Mach number.
     * \return Subsonic split contribution.
     */
    scalar M2Function
    (
        const label side,
        const scalar M
    ) const;

    /**
     * \brief Full AUSM Mach splitting function.
     * \param side Select plus/minus polynomial branch.
     * \param beta Variant-dependent Mach-split coefficient.
     * \param M    Local Mach number.
     * \return Split Mach contribution from one side of the interface.
     */
    scalar M4Function
    (
        const label side,
        const scalar beta,
        const scalar M
    ) const;

    /**
     * \brief AUSM pressure splitting function.
     * \param side Select plus/minus polynomial branch.
     * \param alpha Variant-dependent pressure-split coefficient.
     * \param M     Local Mach number.
     * \return Split pressure weight from one side of the interface.
     */
    scalar P5Function
    (
        const label side,
        const scalar alpha,
        const scalar M
    ) const;

    /**
     * \brief Apply the coefficient set associated with the selected variant.
     * \param variant Runtime-selected AUSM variant.
     */
    void setVariant(const variantType variant);

    /**
     * \brief Compute AUSM mass, momentum, and energy fluxes at one face point.
     *
     * The states are rotated into the face-local frame, the interface Mach
     * number and pressure are reconstructed using the configured AUSM variant,
     * and the momentum flux is rotated back to Cartesian coordinates.
     */
    void calcAUSMFlux
    (
        const scalar rhoL,
        const scalar rhoR,
        const vector& ULv,
        const vector& URv,
        const scalar EL,
        const scalar ER,
        const scalar pL,
        const scalar pR,
        const scalar aL,
        const scalar aR,
        const vector& n,
        scalar& massFlux,
        vector& momentumFlux,
        scalar& energyFlux
    ) const;

public:

    TypeName("AUSM");

    /**
     * \brief Construct from name, dictionary, and DG mesh.
     *
     * The constructor reads the optional \c AUSMCoeffs block and configures
     * the selected AUSM variant.
     */
    dgAUSMFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    virtual ~dgAUSMFluxSolver() = default;

    /**
     * \brief Compute scalar-valued numerical fluxes.
     *
     * This overload is used for mass and energy transport equations and
     * returns the scalar interface flux projected along the face normal.
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;

    /**
     * \brief Compute vector-valued numerical fluxes for momentum transport.
     *
     * The returned quantity is the full Cartesian momentum flux vector at
     * each face Gauss point.
     */
    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;

    /** AUSM stores no cache between calls. */
    virtual void reset() override {}

    /**
     * \brief Read solver coefficients from the dictionary.
     *
     * Supported syntax:
     * \verbatim
     * fluxSolversCoeffs
     * {
     *     AUSMCoeffs
     *     {
     *         variant AUSM3;
     *         Mco     0.01;
     *     }
     * }
     * \endverbatim
     *
     * The \c Mco entry is only used by the \c AUSM3 variant.
     */
    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
