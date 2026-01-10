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
    Foam::idealGas

Description
    Ideal-gas equation of state.

SourceFiles
    idealGas.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_idealGas_H
#define Foam_idealGas_H

#include "eqnOfState.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class idealGas Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::idealGas
 * \brief Ideal-gas equation of state
 *
 * \details
 * idealGas implements the classical ideal-gas equation of state:
 *
 * \f[
 *   p = \rho R T
 * \f]
 *
 * where:
 * - \f$p\f$ is the thermodynamic pressure
 * - \f$\rho\f$ is the density
 * - \f$T\f$ is the temperature
 * - \f$R\f$ is the specific gas constant
 *
 * This EOS:
 * - Assumes calorically or thermally perfect gas behaviour
 * - Is compatible with perfect-gas and thermal-perfect-gas thermo laws
 * - Provides all fundamental thermodynamic relations required by
 *   dgThermoConservative
 *
 * \par Supported relations
 * \li \f$\rho = \frac{p}{R T}\f$
 * \li \f$p = \rho R T\f$
 * \li \f$T = \frac{p}{\rho R}\f$
 *
 * \ingroup grpDGThermo
 */
class idealGas
:
    public eqnOfState
{
private:

    /** Specific gas constant \f$R\f$ [J/(kg·K)] */
    scalar R_;

public:

    //- Runtime type information
    TypeName("idealGas");

    //- Runtime selection table
    declareRunTimeSelectionTable
    (
        autoPtr, idealGas, dictionary,
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh
        ),
        (name, dict, mesh)
    );

    /**
     * \brief Construct from dictionary
     *
     * Reads the gas constant \f$R\f$ from the coefficient dictionary.
     *
     * \param name Model name
     * \param dict Coefficient dictionary
     * \param mesh DG geometric mesh
     */
    idealGas
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /** Virtual destructor */
    virtual ~idealGas() = default;

    /**
     * \brief Return specific gas constant
     *
     * \return Gas constant \f$R\f$
     */
    inline scalar R() const override
    {
        return R_;
    }

    //==============================================================//
    //        Ideal-gas thermodynamic relations
    //==============================================================//

    /**
     * \brief Compute density from pressure and temperature
     *
     * Computes:
     * \f[
     *   \rho = \frac{p}{R T}
     * \f]
     */
    virtual void calcRhoFromPT
    (
        const label cellI,
        const GaussField<scalar>& p,
        const GaussField<scalar>& T,
        GaussField<scalar>& rho
    ) const override;

    /**
     * \brief Compute pressure from density and temperature
     *
     * Computes:
     * \f[
     *   p = \rho R T
     * \f]
     */
    virtual void calcPFromRhoT
    (
        const label cellI,
        const GaussField<scalar>& rho,
        const GaussField<scalar>& T,
        GaussField<scalar>& p
    ) const override;

    /**
     * \brief Compute temperature from pressure and density
     *
     * Computes:
     * \f[
     *   T = \frac{p}{\rho R}
     * \f]
     */
    virtual void calcTFromPRho
    (
        const label cellI,
        const GaussField<scalar>& p,
        const GaussField<scalar>& rho,
        GaussField<scalar>& T
    ) const override;

    //==============================================================//
    //              Capability flags
    //==============================================================//

    /** Ideal gas supports \f$\rho(p,T)\f$ */
    virtual bool canCalcRhoFromPT() const override { return true; }

    /** Ideal gas supports \f$p(\rho,T)\f$ */
    virtual bool canCalcPFromRhoT() const override { return true; }

    /** Ideal gas supports \f$T(p,\rho)\f$ */
    virtual bool canCalcTFromPRho() const override { return true; }

    /** Identify ideal-gas EOS */
    virtual bool isIdealGas() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
