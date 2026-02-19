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
    Foam::transportLaw

Description
    Base class for molecular transport property models.

SourceFiles
    transportLaw.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_transportLaw_H
#define Foam_transportLaw_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "word.H"
#include "GaussField.H"
#include "dgGeomMesh.H"
#include "dgMath.H"
#include "thermoLaw.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class transportLaw Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::transportLaw
 * \brief Base class for DG molecular transport models
 *
 * \details
 * transportLaw defines the interface for evaluating molecular transport
 * properties at DG Gauss points.
 *
 * The class provides Gauss-point-resolved evaluation of:
 * - Dynamic viscosity \f$\mu\f$
 * - Thermal conductivity \f$\kappa\f$
 * - Prandtl number \f$\mathrm{Pr}\f$
 *
 * Transport properties may depend on temperature only,
 * or on both temperature and density, depending on the
 * selected model.
 *
 * Supported transport model categories include:
 * - Constant transport
 * - Sutherland-type ideal-gas transport
 * - Polynomial or power-law transport
 * - Real-gas transport
 * - Kinetic / rarefied-gas transport (VHS/VSS/BGK)
 *
 * transportLaw instances are created via runtime selection
 * and are tightly coupled with a thermoLaw model.
 */
class transportLaw
{
protected:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Model instance name */
        word name_;

        /** Coefficient dictionary (immutable snapshot) */
        const dictionary coeff_;

        /** Reference to DG geometric mesh (non-owning) */
        const dgGeomMesh& mesh_;

        /** Reference to thermodynamic law (non-owning) */
        const thermoLaw& thermo_;

public:

    //====================================================================//
    //                     Runtime type information
    //====================================================================//

    TypeName("transportLaw");

    declareRunTimeSelectionTable
    (
        autoPtr,
        transportLaw,
        dictionary,
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh,
            const thermoLaw& thermo
        ),
        (name, dict, mesh, thermo)
    );

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /**
         * \brief Construct from name, dictionary, mesh, and thermo model
         *
         * \param name Model instance name
         * \param dict Transport-law coefficient dictionary
         * \param mesh DG geometric mesh
         * \param thermo Thermodynamic law
         */
        transportLaw
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh,
            const thermoLaw& thermo
        );

        /** Disallow default construction */
        transportLaw() = delete;

        /** Virtual destructor */
        virtual ~transportLaw() = default;

        /**
         * \brief Construct via runtime selection
         *
         * \param name Model instance name
         * \param dict Transport-law coefficient dictionary
         * \param mesh DG geometric mesh
         * \param thermo Thermodynamic law
         */
        static autoPtr<transportLaw> New
        (
            const word& name,
            const dictionary& dict,
            const dgGeomMesh& mesh,
            const thermoLaw& thermo
        );

    //====================================================================//
    //                         Transport evaluation
    //====================================================================//

        /**
         * \brief Compute dynamic viscosity
         *
         * \param cellI DG cell index
         * \param T Temperature field at Gauss points
         * \param mu Dynamic viscosity field at Gauss points
         */
        virtual void calcMu
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& mu
        ) const = 0;

        /**
         * \brief Compute thermal conductivity
         *
         * \param cellI DG cell index
         * \param T Temperature field at Gauss points
         * \param kappa Thermal conductivity field at Gauss points
         */
        virtual void calcKappa
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& kappa
        ) const = 0;

        /**
         * \brief Compute Prandtl number
         *
         * \param cellI DG cell index
         * \param T Temperature field at Gauss points
         * \param Pr Prandtl number field at Gauss points
         */
        virtual void calcPr
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& Pr
        ) const = 0;

    //====================================================================//
    //                               I/O
    //====================================================================//

        /**
         * \brief Read model coefficients
         *
         * Default implementation performs no action.
         */
        virtual void read();

    //====================================================================//
    //                     Transport-law classification
    //====================================================================//

        /** Return true if constant-transport model */
        virtual bool isConstantTransportLaw() const { return false; }

        /** Return true if Sutherland-type transport model */
        virtual bool isSutherlandTransportLaw() const { return false; }

        /** Return true if polynomial or power-law transport model */
        virtual bool isPolynomialTransportLaw() const { return false; }

        /** Return true if real-gas transport model */
        virtual bool isRealGasTransportLaw() const { return false; }

        /** Return true if kinetic / rarefied-gas transport model */
        virtual bool isKineticTransportLaw() const { return false; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
