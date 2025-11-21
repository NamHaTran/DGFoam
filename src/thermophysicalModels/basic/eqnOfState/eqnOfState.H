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
    Foam::eqnOfState

Description
    Base class for Equation of State (single-phase).
    - Constructible from a coefficients dictionary.
    - Provides rho(T,p) and p(T,rho).

    EOS Classification Notes
    ------------------------

    To ensure compatibility between eqnOfState models and thermoLaw models,
    each EOS provides a set of classification flags. These flags allow higher-
    level thermo models (e.g., constantCp, perfectGasThermo, vibrationalThermo)
    to verify whether the selected EOS is physically compatible with their
    assumptions.

    EOS are grouped into the following categories:

    1. Ideal-gas–type EOS (isIdealGas() = true)
       - EOS of the form p = rho * R * T.
       - Includes classic ideal gas and more advanced thermal-perfect gas
         variants that keep the ideal-gas pressure relation but have Cp(T),
         Cv(T), e(T), h(T) defined elsewhere via thermoLaw.
       - Compatible with most thermoLaw models (constantCp, perfectGasThermo,
         JANAF, vibrational models, polyatomic gas, etc.).

    2. Incompressible EOS (isIncompressible() = true)
       - Density is constant or independent of p and T.
       - Pressure does not participate in thermodynamic relations.
       - Only compatible with thermo models specifically designed for liquids
         or incompressible fluids.

    3. Real-gas EOS (isRealGas() = true)
       - Non-ideal equations of state such as Peng–Robinson, Redlich–Kwong,
         stiffened gas, or Van der Waals.
       - Require thermo models that are derived for real-gas properties.

    4. Kinetic-theory EOS (isKineticModel() = true)
       - EOS derived from kinetic-theory (Boltzmann/ES-BGK), where pressure
         may originate from the stress tensor rather than p = rho*R*T.
       - Used for rarefied flows, micro-scale gas dynamics, or high-Kn models.

    By default, all flags return false in the base class. Individual EOS
    implementations override the appropriate flag(s) to identify their group.

    Thermo models should check these flags in their constructors to ensure
    that the selected EOS is physically consistent with the thermo model.

SourceFiles
    eqnOfState.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_eqnOfState_H
#define Foam_eqnOfState_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "word.H"
#include "dgField.H"
#include "dgGeomMesh.H"
#include "dgMath.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class eqnOfState Declaration
\*---------------------------------------------------------------------------*/

class eqnOfState
{
protected:

    //- Model name (for logging/debug)
    word name_;

    //- Coefficients dictionary (owned/const view)
    const dictionary dict_;

    //- Reference to DG mesh (non-owning)
    const dgGeomMesh& mesh_;

public:

    //- Runtime type information
    TypeName("eqnOfState");

    //- Runtime selection: construct from dictionary
    declareRunTimeSelectionTable
    (
        autoPtr, eqnOfState, dictionary,
        (const word& name, const dictionary& dict, const dgGeomMesh& mesh),
        (name, dict, mesh)
    );

    //- Construct from model name and coefficients dictionary
    eqnOfState
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //- Destructor
    virtual ~eqnOfState() = default;

    //- Factory: create eqnOfState from dictionary
    static autoPtr<eqnOfState> New
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    //====================================================================//
    //        General EOS API (all EOS classes must provide these)
    //====================================================================//

    //- Return gas constant R (returns 0 for incompressible models)
    virtual scalar R() const = 0;


    //====================================================================//
    //   Fundamental thermodynamic relations (general Gaussian operations)
    //
    //   These functions receive GaussField for p, T, rho and return the
    //   corresponding Gaussian field for the selected cell.
    //
    //   NOTE:
    //     Models that do not support a relation simply set the flag below
    //     and dgThermo will not call these functions.
    //====================================================================//

    //--- Compute rho = rho(p, T)
    virtual void calcRhoFromPT
    (
        const label cellI,
        const GaussField<scalar>& p,
        const GaussField<scalar>& T,
        GaussField<scalar>& rho
    ) const
    {
        FatalErrorInFunction
            << "calcRhoFromPT() not supported by this EOS model"
            << exit(FatalError);
    }

    //--- Compute p = p(rho, T)
    virtual void calcPFromRhoT
    (
        const label cellI,
        const GaussField<scalar>& rho,
        const GaussField<scalar>& T,
        GaussField<scalar>& p
    ) const
    {
        FatalErrorInFunction
            << "calcPFromRhoT() not supported by this EOS model"
            << exit(FatalError);
    }

    //--- Compute T = T(p, rho)
    virtual void calcTFromPRho
    (
        const label cellI,
        const GaussField<scalar>& p,
        const GaussField<scalar>& rho,
        GaussField<scalar>& T
    ) const
    {
        FatalErrorInFunction
            << "calcTFromPRho() not supported by this EOS model"
            << exit(FatalError);
    }

    //====================================================================//
    //                EOS capability flags (very important)
    //
    //   Allow dgThermo to detect which relation is supported.
    //   PerfectGas:   all true
    //   Incompressible: all false, except rhoConst()
    //====================================================================//

    virtual bool canCalcRhoFromPT()  const { return false; }
    virtual bool canCalcPFromRhoT()  const { return false; }
    virtual bool canCalcTFromPRho()  const { return false; }

    // Equation of state groups

        //- Ideal-gas–based EOS
        virtual bool isIdealGas() const { return false; }

        //- Incompressible EOS
        virtual bool isIncompressible() const { return false; }

        //- Real-gas EOS
        virtual bool isRealGas() const { return false; }

        //- Special EOS (kinetic/ES-BGK/stress-based)
        virtual bool isThermalPerfectGas() const { return false; }

};

} // End namespace Foam

#endif

// ************************************************************************* //
