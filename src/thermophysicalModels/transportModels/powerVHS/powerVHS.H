/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Foam::powerVHS

Description
    Variable-Hard-Sphere (power-law) transport model.

SourceFiles
    powerVHS.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_powerVHS_H
#define Foam_powerVHS_H

#include "transportLaw.H"
#include "word.H"
#include "scalar.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class powerVHS Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::powerVHS
 * \brief Variable-Hard-Sphere (VHS) power-law viscosity model
 *
 * \details
 * powerVHS implements the Variable-Hard-Sphere (VHS) model with a power-law
 * temperature dependence of viscosity:
 *
 * \f[
 *   \mu(T) = \mu_{\mathrm{ref}}
 *            \left(\frac{T}{T_{\mathrm{ref}}}\right)^{\omega}
 * \f]
 *
 * The reference viscosity is computed from kinetic theory as:
 *
 * \f[
 *   \mu_{\mathrm{ref}} =
 *   \frac{15 \sqrt{\pi m k_B T_{\mathrm{ref}}}}
 *        {2 \pi d_{\mathrm{ref}}^2
 *         (5 - 2\omega)(7 - 2\omega)}
 * \f]
 *
 * where:
 * \li \f$m\f$ is the molecular mass
 * \li \f$k_B\f$ is the Boltzmann constant
 * \li \f$d_{\mathrm{ref}}\f$ is the reference molecular diameter
 * \li \f$T_{\mathrm{ref}}\f$ is the reference temperature
 * \li \f$\omega\f$ is the VHS power-law exponent
 *
 * The Prandtl number is assumed constant:
 * \f[
 *   \mathrm{Pr}(T) = \mathrm{Pr}_0
 * \f]
 *
 * Thermal conductivity is computed in the base transportLaw as:
 * \f[
 *   \kappa(T) = \frac{\mu(T)\,C_p(T)}{\mathrm{Pr}(T)}
 * \f]
 *
 * \par Notes
 * \li Molecular mass is read in \f$\mathrm{g/mol}\f$ and converted internally
 *     to \f$\mathrm{kg/molecule}\f$
 * \li \f$d_{\mathrm{ref}}\f$ is specified in \f$\mathrm{m}\f$
 * \li \f$T_{\mathrm{ref}}\f$ is specified in \f$\mathrm{K}\f$
 * \li \f$k_B\f$ is specified in \f$\mathrm{J/K}\f$
 * \li \f$N_A\f$ is specified in \f$\mathrm{mol^{-1}}\f$
 *
 * \par Example
 * \code
 * transportCoeffs
 * {
 *     type     powerVHS;
 *
 *     // Required
 *     molMass  28.0134;     // [g/mol] N2
 *     dRef     3.7e-10;     // [m]
 *     TRef     300;         // [K]
 *     omega    0.74;        // [-]
 *     Pr       0.72;        // [-]
 *
 *     // Optional overrides
 *     // kB    1.380649e-23;   // [J/K]
 *     // NA    6.02214076e23;  // [1/mol]
 * }
 * \endcode
 *
 * \ingroup grpDGThermo
 */
class powerVHS
:
    public transportLaw
{
private:

    // Coefficients read from dictionary

        /** Molecular mass \f$m\f$ [g/mol] (converted internally) */
        scalar molMass_;

        /** Reference molecular diameter \f$d_{\mathrm{ref}}\f$ [m] */
        scalar dRef_;

        /** Reference temperature \f$T_{\mathrm{ref}}\f$ [K] */
        scalar TRef_;

        /** VHS power-law exponent \f$\omega\f$ [-] */
        scalar omega_;

        /** Constant Prandtl number \f$\mathrm{Pr}_0\f$ [-] */
        scalar Pr0_;

        /** Boltzmann constant \f$k_B\f$ [J/K] (optional override) */
        scalar kB_;

        /** Avogadro constant \f$N_A\f$ [mol^{-1}] (optional override) */
        scalar NA_;

    // Derived quantities

        /** Reference viscosity \f$\mu_{\mathrm{ref}}\f$ at \f$T_{\mathrm{ref}}\f$
            [kg/(m·s)] */
        scalar muRef_;

public:

    //- Runtime type information
    TypeName("powerVHS");

    /**
     * \brief Construct from dictionary
     *
     * Reads transport coefficients and pre-computes
     * the reference viscosity \f$\mu_{\mathrm{ref}}\f$.
     */
    powerVHS
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    /** Virtual destructor */
    virtual ~powerVHS() override = default;


    //====================================================================//
    //                   Transport-law interface
    //====================================================================//

    /**
     * \brief Compute dynamic viscosity
     *
     * Computes:
     * \f[
     *   \mu(T) = \mu_{\mathrm{ref}}
     *            \left(\frac{T}{T_{\mathrm{ref}}}\right)^{\omega}
     * \f]
     */
    virtual void calcMu
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& mu
    ) const override;

    /**
     * \brief Compute thermal conductivity
     *
     * Uses:
     * \f[
     *   \kappa(T) = \frac{\mu(T)\,C_p(T)}{\mathrm{Pr}_0}
     * \f]
     */
    virtual void calcKappa
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& kappa
    ) const override;

    /**
     * \brief Return Prandtl number
     *
     * Returns constant:
     * \f[
     *   \mathrm{Pr}(T) = \mathrm{Pr}_0
     * \f]
     */
    virtual void calcPr
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& Pr
    ) const override;

    /**
     * \brief Read coefficients and pre-compute reference viscosity
     */
    virtual void read() override;

    /**
     * \brief Identify kinetic-theory transport model
     */
    virtual bool isKineticTransportLaw() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
