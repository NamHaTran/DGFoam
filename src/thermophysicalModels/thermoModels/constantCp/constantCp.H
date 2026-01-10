/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::constantCp

Description
    Thermodynamic model with constant specific heat at constant pressure.

SourceFiles
    constantCp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_constantCp_H
#define Foam_constantCp_H

#include "thermoLaw.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class constantCp Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::constantCp
 * \brief Perfect-gas thermodynamic model with constant \f$C_p\f$
 *
 * \details
 * constantCp implements a thermodynamic model in which the specific heat
 * at constant pressure is assumed to be constant.
 *
 * The model reads:
 * \li \f$C_p\f$ from the \c thermodynamics sub-dictionary
 * \li Molecular weight from the \c specie sub-dictionary
 *
 * The following thermodynamic properties are then derived:
 * \li Gas constant: \f$R = \frac{R_u}{W}\f$
 * \li Specific heat at constant volume: \f$C_v = C_p - R\f$
 * \li Heat-capacity ratio: \f$\gamma = \frac{C_p}{C_v}\f$
 *
 * Energy relations:
 * \li Internal energy: \f$e(T) = C_v T\f$
 * \li Enthalpy: \f$h(T) = C_p T\f$
 *
 * Speed of sound:
 * \f[
 *   a = \sqrt{\gamma R T}
 * \f]
 *
 * \par Applicability
 * \li Ideal-gas equation of state
 * \li Perfect-gas thermodynamics
 * \li Compressible DG solvers
 *
 * \par Limitations
 * \li Does not include formation enthalpy
 * \li Not suitable for real-gas or reacting-flow thermodynamics
 *
 * \ingroup grpDGThermo
 */
class constantCp
:
    public thermoLaw
{
private:

    /** Constant specific heat at constant pressure \f$C_p\f$ [J/(kg·K)] */
    scalar Cp_;

    /** Specific gas constant \f$R\f$ [J/(kg·K)] */
    scalar R_;

    /** Constant specific heat at constant volume \f$C_v\f$ [J/(kg·K)] */
    scalar Cv_;

    /** Heat-capacity ratio \f$\gamma = C_p / C_v\f$ */
    scalar gamma_;

public:

    //- Runtime type information
    TypeName("constantCp");

    /**
     * \brief Construct from root dictionary
     *
     * Reads:
     * \li \c Cp from \c thermodynamics sub-dictionary
     * \li Molecular weight from \c specie sub-dictionary
     *
     * \param name     Model name
     * \param rootDict Root dictionary
     * \param mesh     DG geometric mesh
     * \param eos      Equation-of-state model
     */
    constantCp
    (
        const word& name,
        const dictionary& rootDict,
        const dgGeomMesh& mesh,
        const eqnOfState& eos
    );

    /** Virtual destructor */
    virtual ~constantCp() = default;

    //====================================================================//
    //                   Thermodynamic interface
    //====================================================================//

    /**
     * \brief Compute heat capacity at constant pressure
     *
     * Returns constant \f$C_p\f$.
     */
    virtual void calcCp
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& Cp
    ) const override;

    /**
     * \brief Compute heat capacity at constant volume
     *
     * Returns constant \f$C_v\f$.
     */
    virtual void calcCv
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& Cv
    ) const override;

    /**
     * \brief Compute heat-capacity ratio
     *
     * Computes:
     * \f[
     *   \gamma = \frac{C_p}{C_v}
     * \f]
     */
    virtual void calcGamma
    (
        const label cellI,
        const GaussField<scalar>& CpField,
        const GaussField<scalar>& CvField,
        GaussField<scalar>& gamma
    ) const override;

    /**
     * \brief Compute temperature from energy
     *
     * Inverts:
     * \li \f$e = C_v T\f$ (internal-energy formulation)
     * \li \f$h = C_p T\f$ (enthalpy formulation)
     */
    virtual void calcT
    (
        const label cellI,
        const GaussField<scalar>& energy,
        GaussField<scalar>& T
    ) const override;

    /**
     * \brief Compute internal energy
     *
     * Computes:
     * \f[
     *   e(T) = C_v T
     * \f]
     */
    virtual void calcInternalE
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& e
    ) const override;

    /**
     * \brief Compute enthalpy
     *
     * Computes:
     * \f[
     *   h(T) = C_p T
     * \f]
     */
    virtual void calcH
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& h
    ) const override;

    /**
     * \brief Compute speed of sound
     *
     * Computes:
     * \f[
     *   a = \sqrt{\gamma R T}
     * \f]
     */
    virtual void calcSpeedOfSound
    (
        const label cellI,
        const GaussField<scalar>& T,
        const GaussField<scalar>& gamma,
        GaussField<scalar>& a
    ) const override;

    /**
     * \brief Identify perfect-gas thermodynamic model
     */
    virtual bool isPerfectGasThermo() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
