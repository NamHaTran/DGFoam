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
    Foam::pengRobinson

Description
    Peng-Robinson equation of state for real-gas thermodynamic closure.

SourceFiles
    pengRobinson.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pengRobinson_H
#define Foam_pengRobinson_H

#include "eqnOfState.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class pengRobinson Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::pengRobinson
 * \brief Peng-Robinson real-gas equation of state for single-species flows
 *
 * \par Declaration syntax
 * \code
 * class pengRobinson
 * :
 *     public eqnOfState
 * {
 *     // ...
 * };
 * \endcode
 *
 * \details
 * This class implements the Peng-Robinson cubic equation of state (EOS)
 * for compressible single-phase fluids:
 *
 * \f[
 *   p = \frac{RT}{1/\rho - b}
 *     - \frac{a \alpha(T)}{1/\rho^2 + 2b/\rho - b^2}
 * \f]
 *
 * with the standard Peng-Robinson coefficients:
 *
 * \f[
 *   a = 0.45724 \frac{R^2 T_c^2}{P_c},
 *   \qquad
 *   b = 0.07780 \frac{R T_c}{P_c}
 * \f]
 *
 * and temperature correction:
 *
 * \f[
 *   \alpha(T)
 *   =
 *   \left[1 + f_\omega
 *   \left(1 - \sqrt{T/T_c}\right)\right]^2
 * \f]
 *
 * where:
 *
 * \f[
 *   f_\omega = 0.37464 + 1.54226\omega - 0.26992\omega^2
 * \f]
 *
 * The implementation exposes the conservative real-gas closures needed by
 * \c rhoBasedConservative:
 *
 * - \f$\rho = \rho(p,T)\f$
 * - \f$p = p(\rho,T)\f$
 * - \f$T = T(p,\rho)\f$
 * - \f$T = T(\rho,e)\f$
 * - \f$p = p(\rho,e)\f$
 * - \f$e = e(\rho,p)\f$
 * - \f$e = e(\rho,T)\f$
 * - \f$a = a(\rho,e)\f$
 *
 * \par Thermodynamic assumptions
 * The Peng-Robinson EOS here is used together with a reference caloric model
 * from \c thermoLaw, currently expected to be \c constantCp in DGFoam's
 * real-gas conservative workflow. The EOS itself owns the density/pressure/
 * temperature non-ideal closure and reconstructs the internal energy and
 * speed of sound required by the solver.
 *
 * \par Required dictionary entries
 * The model reads data from the usual \c mixture dictionary:
 *
 * \code
 * mixture
 * {
 *     specie
 *     {
 *         molWeight       44.01;
 *     }
 *
 *     thermodynamics
 *     {
 *         Cp              846.0;
 *     }
 *
 *     equationOfState
 *     {
 *         Tcrit           304.1282;    // or Tc
 *         Pcrit           7.3773e6;    // or Pc
 *         AcentricFactor  0.22394;     // or omega
 *     }
 * }
 * \endcode
 *
 * \ingroup grpDGThermo
 */
class pengRobinson
:
    public eqnOfState
{
private:

    /** Specific gas constant \f$R\f$ [J/(kg.K)] */
    scalar R_;

    /** Critical temperature \f$T_c\f$ [K] */
    scalar Tc_;

    /** Critical pressure \f$P_c\f$ [Pa] */
    scalar Pc_;

    /** Pitzer acentric factor \f$\omega\f$ [-] */
    scalar omega_;

    /** Peng-Robinson temperature-shape coefficient \f$f_\omega\f$ [-] */
    scalar fw_;

    /** Attraction coefficient \f$a\f$ in the cubic EOS */
    scalar aCoeff_;

    /** Co-volume coefficient \f$b\f$ in the cubic EOS */
    scalar bCoeff_;

    /** Reference constant-pressure heat capacity used for caloric closure */
    scalar cpIdeal_;

    /** Reference constant-volume heat capacity used for caloric closure */
    scalar cvIdeal_;

    /** Reference ideal-gas heat-capacity ratio */
    scalar gammaIdeal_;

    /**
     * \brief Evaluate Peng-Robinson temperature correction \f$\alpha(T)\f$
     *
     * \param T Temperature [K]
     * \return \f$\alpha(T)\f$
     */
    scalar alpha(const scalar T) const;

    /**
     * \brief Evaluate the logarithmic density term used in real-gas energy
     * corrections
     *
     * \param rho Density [kg/m^3]
     * \return Logarithmic term appearing in the Peng-Robinson departure
     *         functions
     */
    scalar logTerm(const scalar rho) const;

    /**
     * \brief Recover temperature from density and pressure
     *
     * \details
     * Solves the quadratic form in \f$\sqrt{T}\f$ that results from the
     * Peng-Robinson EOS when \f$\rho\f$ and \f$p\f$ are prescribed.
     *
     * \param rho Density [kg/m^3]
     * \param p   Pressure [Pa]
     * \return Temperature [K]
     */
    scalar calcTemperatureFromRhoP(const scalar rho, const scalar p) const;

    /**
     * \brief Evaluate \f$\partial p / \partial e \vert_\rho\f$
     *
     * \param rho Density [kg/m^3]
     * \param e   Specific internal energy [J/kg]
     * \return Pressure derivative with respect to internal energy at
     *         constant density
     */
    scalar calcDPDeRho(const scalar rho, const scalar e) const;

    /**
     * \brief Evaluate \f$\partial p / \partial \rho \vert_e\f$
     *
     * \param rho Density [kg/m^3]
     * \param e   Specific internal energy [J/kg]
     * \return Pressure derivative with respect to density at constant
     *         internal energy
     */
    scalar calcDPDrhoE(const scalar rho, const scalar e) const;

public:

    TypeName("pengRobinson");

    using eqnOfState::calcRhoFromPT;
    using eqnOfState::calcPFromRhoT;
    using eqnOfState::calcTFromPRho;

    /**
     * \brief Construct from EOS name, coefficients dictionary, and DG mesh
     *
     * \details
     * The constructor reads:
     * - \c specie/molWeight
     * - \c thermodynamics/Cp
     * - \c equationOfState/Tcrit (or \c Tc)
     * - \c equationOfState/Pcrit (or \c Pc)
     * - \c equationOfState/AcentricFactor (or \c omega)
     *
     * It then computes the Peng-Robinson constants \f$a\f$, \f$b\f$, and
     * \f$f_\omega\f$, along with the reference ideal-gas caloric constants
     * used in the conservative real-gas closure.
     *
     * \param name Model name from the runtime-selection table
     * \param dict Root mixture dictionary
     * \param mesh DG geometric mesh
     */
    pengRobinson
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    virtual ~pengRobinson() = default;

    inline scalar R() const override
    {
        return R_;
    }

    /**
     * \brief Compute density from pressure and temperature
     *
     * \details
     * Solves the Peng-Robinson compressibility equation for \f$Z\f$ with a
     * Newton iteration and returns \f$\rho = p/(ZRT)\f$.
     *
     * \param p Pressure [Pa]
     * \param T Temperature [K]
     * \return Density [kg/m^3]
     */
    virtual scalar calcRhoFromPT
    (
        const scalar p,
        const scalar T
    ) const override;

    /**
     * \brief Compute pressure from density and temperature
     *
     * \param rho Density [kg/m^3]
     * \param T   Temperature [K]
     * \return Pressure [Pa]
     */
    virtual scalar calcPFromRhoT
    (
        const scalar rho,
        const scalar T
    ) const override;

    /**
     * \brief Compute temperature from pressure and density
     *
     * \param p   Pressure [Pa]
     * \param rho Density [kg/m^3]
     * \return Temperature [K]
     */
    virtual scalar calcTFromPRho
    (
        const scalar p,
        const scalar rho
    ) const override;

    /**
     * \brief Compute temperature from density and internal energy
     *
     * \details
     * Uses the Peng-Robinson real-gas internal-energy relation and solves the
     * resulting quadratic equation in \f$\sqrt{T}\f$.
     *
     * \param rho Density [kg/m^3]
     * \param e   Specific internal energy [J/kg]
     * \return Temperature [K]
     */
    virtual scalar calcTFromRhoE
    (
        const scalar rho,
        const scalar e
    ) const override;

    /**
     * \brief Compute pressure from density and internal energy
     *
     * \param rho Density [kg/m^3]
     * \param e   Specific internal energy [J/kg]
     * \return Pressure [Pa]
     */
    virtual scalar calcPFromRhoE
    (
        const scalar rho,
        const scalar e
    ) const override;

    /**
     * \brief Compute internal energy from density and pressure
     *
     * \param rho Density [kg/m^3]
     * \param p   Pressure [Pa]
     * \return Specific internal energy [J/kg]
     */
    virtual scalar calcEFromRhoP
    (
        const scalar rho,
        const scalar p
    ) const override;

    /**
     * \brief Compute internal energy from density and temperature
     *
     * \param rho Density [kg/m^3]
     * \param T   Temperature [K]
     * \return Specific internal energy [J/kg]
     */
    virtual scalar calcEFromRhoT
    (
        const scalar rho,
        const scalar T
    ) const override;

    /**
     * \brief Compute speed of sound from density and internal energy
     *
     * \details
     * Evaluates the thermodynamic sound speed from the identity
     * based on \f$\partial p/\partial \rho \vert_e\f$,
     * \f$\partial p/\partial e \vert_\rho\f$, and enthalpy.
     *
     * \param rho Density [kg/m^3]
     * \param e   Specific internal energy [J/kg]
     * \return Speed of sound [m/s]
     */
    virtual scalar calcAFromRhoE
    (
        const scalar rho,
        const scalar e
    ) const override;

    /**
     * \brief Compute \f$\nabla p\f$ for the explicit PR relation p(rho,T).
     */
    virtual vector calcGradPFromRhoT
    (
        const scalar rho,
        const scalar T,
        const vector& gradRho,
        const vector& gradT
    ) const override;

    virtual bool canCalcRhoFromPT() const override { return true; }
    virtual bool canCalcPFromRhoT() const override { return true; }
    virtual bool canCalcTFromPRho() const override { return true; }
    virtual bool canCalcTFromRhoE() const override { return true; }
    virtual bool canCalcPFromRhoE() const override { return true; }
    virtual bool canCalcEFromRhoP() const override { return true; }
    virtual bool canCalcEFromRhoT() const override { return true; }
    virtual bool canCalcAFromRhoE() const override { return true; }

    virtual bool isRealGas() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
