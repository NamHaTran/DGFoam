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
    Foam::sensibleInternalEnergy

Description
    Sensible internal-energy model for DG thermodynamics.

SourceFiles
    sensibleInternalEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sensibleInternalEnergy_H
#define Foam_sensibleInternalEnergy_H

#include "energy.H"
#include "word.H"
#include "dictionary.H"
#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class sensibleInternalEnergy Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::sensibleInternalEnergy
 * \brief Sensible internal-energy formulation \f$e = C_v T\f$
 *
 * \details
 * sensibleInternalEnergy implements an energy model in which the conserved
 * thermodynamic variable is the sensible internal energy.
 *
 * The internal energy is defined as:
 * \f[
 *   e(T) = C_v(T)\,T
 * \f]
 *
 * Enthalpy is recovered consistently as:
 * \f[
 *   h(T) = e(T) + \frac{p}{\rho}
 * \f]
 *
 * This model:
 * - Evolves sensible internal energy only
 * - Excludes chemical formation energy
 * - Is suitable for compressible DG solvers
 * - Is compatible with perfect-gas and thermal-perfect-gas thermodynamics
 *
 * \par Usage
 * sensibleInternalEnergy is intended for conservative DG formulations
 * where the energy equation is written in terms of internal energy.
 *
 * \ingroup grpDGThermo
 */
class sensibleInternalEnergy
:
    public energy
{
protected:

    /** Optional coefficient dictionary */
    dictionary coeffDict_;

public:

    //- Runtime type information
    TypeName("sensibleInternalEnergy");

    /**
     * \brief Construct from dictionary
     *
     * \param name   Model name
     * \param dict   Coefficient dictionary
     * \param mesh   DG geometric mesh
     * \param thermo Thermodynamic law
     */
    sensibleInternalEnergy
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    // * * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

    /**
     * \brief Compute enthalpy from temperature
     *
     * Computes sensible enthalpy as:
     * \f[
     *   h(T) = e(T) + \frac{p}{\rho}
     * \f]
     *
     * \param cellI Cell index
     * \param T     Temperature field at Gauss points
     * \param h     Enthalpy field to be written
     */
    virtual void calcEnthalpy
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& h
    ) const override;

    /**
     * \brief Compute internal energy from temperature
     *
     * Computes sensible internal energy:
     * \f[
     *   e(T) = C_v(T)\,T
     * \f]
     *
     * \param cellI Cell index
     * \param T     Temperature field at Gauss points
     * \param e     Internal energy field to be written
     */
    virtual void calcEnergy
    (
        const label cellI,
        const GaussField<scalar>& T,
        GaussField<scalar>& e
    ) const override;

    /**
     * \brief Compute temperature from internal energy
     *
     * Inverts the relation:
     * \f[
     *   e = C_v(T)\,T
     * \f]
     *
     * \param cellI Cell index
     * \param e     Internal energy field at Gauss points
     * \param T     Temperature field to be written
     */
    virtual void calcTfromEnergy
    (
        const label cellI,
        const GaussField<scalar>& e,
        GaussField<scalar>& T
    ) const override;

    /**
     * \brief Identify internal-energy formulation
     *
     * \return True, since this model is internal-energy based
     */
    virtual bool energyInternal() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
