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

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCompressibleRiemannInvariantBoundaryField_H
#define Foam_dgCompressibleRiemannInvariantBoundaryField_H

#include "dgCompressibleBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class dgCompressibleRiemannInvariantBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleRiemannInvariantBoundaryField
 * \brief Far-field compressible DG boundary condition based on Riemann invariants.
 *
 * \details
 * Combines interior and freestream data through a one-dimensional
 * characteristic model in the outward normal direction \f$\mathbf{n}\f$.
 *
 * The implementation first reconstructs the interior primitive state:
 *
 * \f[
 *     \mathbf{U}^- = \frac{\rho\mathbf{U}^-}{\rho^-}, \qquad
 *     V_n^- = \mathbf{n}\cdot \mathbf{U}^-,
 * \qquad
 *     V_{n,\infty} = \mathbf{n}\cdot \mathbf{U}_{\infty},
 * \f]
 *
 * \f[
 *     k^- = \frac{1}{2}\left|\mathbf{U}^-\right|^2, \qquad
 *     h_e^- = \frac{E^-}{\rho^-} - k^-,
 * \qquad
 *     T^- = T(h_e^-),
 * \qquad
 *     p^- = p(\rho^-, T^-).
 * \f]
 *
 * With \f$\gamma=\gamma_{\infty}\f$, the interior and freestream sound speeds are
 *
 * \f[
 *     c^- = \sqrt{\frac{\gamma p^-}{\rho^-}}, \qquad
 *     c_{\infty} = c(T_{\infty}, \gamma),
 * \f]
 *
 * and the Riemann invariants are
 *
 * \f[
 *     R^+ = V_n + \frac{2c}{\gamma-1}, \qquad
 *     R^- = V_n - \frac{2c}{\gamma-1}.
 * \f]
 *
 * The boundary model selects these invariants as follows:
 *
 * \f[
 * \text{subsonic:}\quad
 * R^+ = V_n^- + \frac{2c^-}{\gamma-1}, \qquad
 * R^- = V_{n,\infty} - \frac{2c_{\infty}}{\gamma-1},
 * \f]
 *
 * \f[
 * \text{supersonic inflow:}\quad
 * R^{\pm} = V_{n,\infty} \pm \frac{2c_{\infty}}{\gamma-1},
 * \f]
 *
 * \f[
 * \text{supersonic outflow:}\quad
 * R^{\pm} = V_n^- \pm \frac{2c^-}{\gamma-1}.
 * \f]
 *
 * From these invariants,
 *
 * \f[
 *     V_{n,b} = \frac{1}{2}\left(R^+ + R^-\right), \qquad
 *     c_b = \frac{\gamma-1}{4}\left(R^+ - R^-\right).
 * \f]
 *
 * The tangential velocity and entropy are taken from the freestream for
 * inflow and from the interior state for outflow:
 *
 * \f[
 *     \mathbf{U}_b =
 *     \begin{cases}
 *         \mathbf{U}_{\infty} + (V_{n,b}-V_{n,\infty})\mathbf{n},
 *         & V_n^- \le 0, \\
 *         \mathbf{U}^- + (V_{n,b}-V_n^-)\mathbf{n},
 *         & V_n^- > 0,
 *     \end{cases}
 * \f]
 *
 * \f[
 *     s_b =
 *     \begin{cases}
 *         s_{\infty}, & V_n^- \le 0, \\
 *         \dfrac{p^-}{(\rho^-)^{\gamma}}, & V_n^- > 0.
 *     \end{cases}
 * \f]
 *
 * The final ghost state is recovered from
 *
 * \f[
 *     \rho^+ = \left(\frac{c_b^2}{\gamma s_b}\right)^{\frac{1}{\gamma-1}},
 *     \qquad
 *     p_b = \frac{\rho^+ c_b^2}{\gamma},
 *     \qquad
 *     \rho\mathbf{U}^+ = \rho^+ \mathbf{U}_b,
 * \f]
 *
 * with the total energy assembled consistently with the selected internal-
 * energy or enthalpy formulation in the thermodynamic model.
 */
class dgCompressibleRiemannInvariantBoundaryField
:
    public dgCompressibleBoundaryField
{
private:

    /** \brief Freestream pressure. */
    scalar pInf_;

    /** \brief Freestream temperature. */
    scalar TInf_;

    /** \brief Freestream velocity. */
    vector UInf_;

    /** \brief Freestream density. */
    scalar rhoInf_;

    /** \brief Reference heat-capacity ratio used by the characteristic model. */
    scalar gammaInf_;

    /** \brief Freestream speed of sound. */
    scalar cInf_;

    /** \brief Freestream isentropic invariant \f$s_{\infty}=p_{\infty}/\rho_{\infty}^{\gamma}\f$. */
    scalar sInf_;

public:

    TypeName("compressible::RiemannInvariantBC");

    dgCompressibleRiemannInvariantBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dgThermoConservative& thermo,
        const dictionary& dict
    );

    /** \brief Default destructor. */
    virtual ~dgCompressibleRiemannInvariantBoundaryField() = default;

    /**
     * \brief Compute the ghost state from characteristic matching.
     *
     * \param cellID Owner cell index.
     * \param faceLocalID Local face index on the owner cell.
     * \param localGauss Local Gauss-point index on the face.
     * \param n Outward unit normal of the owner cell.
     * \param rhoMinus Interior density.
     * \param rhoUMinus Interior momentum.
     * \param EMinus Interior total energy.
     * \param rhoPlus Ghost density.
     * \param rhoUPlus Ghost momentum.
     * \param EPlus Ghost total energy.
     */
    virtual void updateGhostState
    (
        const label cellID,
        const label faceLocalID,
        const label localGauss,
        const vector& n,
        const scalar rhoMinus,
        const vector& rhoUMinus,
        const scalar EMinus,
        scalar& rhoPlus,
        vector& rhoUPlus,
        scalar& EPlus
    ) const override;

    /** \brief Check that the boundary condition is applied on a standard patch. */
    virtual void checkPatchType() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
