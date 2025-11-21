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
    Thermodynamic model with constant Cp.
    - Reads Cp from 'thermodynamics' sub-dictionary
    - Reads molWeight from 'specie' sub-dictionary
    - Computes R = Ru / molWeight, Cv = Cp - R, gamma = Cp/Cv
    - Provides cp(T), h(T), e(T) and speed of sound a(T) = sqrt(gamma*R*T)

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

class constantCp
:
    public thermoLaw
{
private:

    //- Constant specific heat at constant pressure [J/(kg.K)]
    scalar Cp_;

    //- Specific gas constant (computed) [J/(kg.K)]
    scalar R_;

    //- Constant specific heat at constant volume (computed) [J/(kg.K)]
    scalar Cv_;

    //- Heat capacity ratio gamma = Cp/Cv [-]
    scalar gamma_;

public:

    //- Runtime type information
    TypeName("constantCp");

    //- Construct from root dictionary (reads subDicts 'thermodynamics' and 'specie')
    constantCp
    (
        const word& name, 
        const dictionary& rootDict, 
        const dgGeomMesh& mesh,
        const eqnOfState& eos
    );

    //- Destructor
    virtual ~constantCp() = default;

    // Member Functions

        //- Cp(T)
        virtual void calcCp
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& Cp
        ) const override;

        //- Cv(T)
        virtual void calcCv
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& Cv
        ) const override;

        //- gamma(T) = Cp / Cv
        virtual void calcGamma
        (
            const label cellI,
            const GaussField<scalar>& CpField,
            const GaussField<scalar>& CvField,
            GaussField<scalar>& gamma
        ) const override;

        //- T from energy (e or h depending on energy model)
        virtual void calcT
        (
            const label cellI,
            const GaussField<scalar>& energy,
            GaussField<scalar>& T
        ) const override;

        //- e(T)
        virtual void calcInternalE
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& e
        ) const override;

        //- h(T)
        virtual void calcH
        (
            const label cellI,
            const GaussField<scalar>& T,
            GaussField<scalar>& h
        ) const override;

        //- Calculate speedOfSound
        virtual void calcSpeedOfSound
        (
            const label cellI,
            const GaussField<scalar>& T,
            const GaussField<scalar>& gamma,
            GaussField<scalar>& a
        ) const override;

        //- Set group
        //- This model is in PerfectGasThermo
        virtual bool isPerfectGasThermo() const override { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
