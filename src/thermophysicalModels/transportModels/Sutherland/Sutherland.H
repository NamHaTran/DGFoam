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
    Foam::Sutherland

Description
    Classical Sutherland-law transport model.

SourceFiles
    Sutherland.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Sutherland_H
#define Foam_Sutherland_H

#include "transportLaw.H"
#include "word.H"
#include "scalar.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Sutherland Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::Sutherland
 * \brief Classical Sutherland-law viscosity model
 *
 * \details
 * Implements the classical Sutherland law for temperature-dependent
 * dynamic viscosity:
 *
 * \f[
 *   \mu(T) = A_s \frac{T^{3/2}}{T + S}
 * \f]
 *
 * The Prandtl number is assumed constant:
 * \f[
 *   \mathrm{Pr}(T) = \mathrm{Pr}_0
 * \f]
 *
 * Thermal conductivity is evaluated as:
 * \f[
 *   \kappa(T) = \frac{\mu(T)\,C_p(T)}{\mathrm{Pr}_0}
 * \f]
 *
 * All evaluations are performed on per-cell GaussField data.
 *
 * \ingroup grpDGThermo
 */
class Sutherland
:
    public transportLaw
{
private:

    // Coefficients read from dictionary

        /** Sutherland coefficient \f$A_s\f$
            [kg/(m·s·\sqrt{K})] */
        scalar As_;

        /** Sutherland temperature \f$S\f$ [K] */
        scalar S_;

        /** Constant Prandtl number \f$\mathrm{Pr}_0\f$ [-] */
        scalar Pr0_;

public:

    //- Runtime type information
    TypeName("Sutherland");

    /**
     * \brief Construct from dictionary
     *
     * Reads Sutherland coefficients and constant Prandtl number.
     */
    Sutherland
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh,
        const thermoLaw& thermo
    );

    /** Virtual destructor */
    virtual ~Sutherland() = default;


    //====================================================================//
    //                   Transport-law interface
    //====================================================================//

    /**
     * \brief Compute dynamic viscosity
     *
     * Evaluates:
     * \f[
     *   \mu(T) = A_s \frac{T^{3/2}}{T + S}
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
     * \brief Read coefficients from dictionary
     */
    virtual void read() override;

    /**
     * \brief Identify Sutherland transport model
     */
    virtual bool isSutherlandTransportLaw() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
