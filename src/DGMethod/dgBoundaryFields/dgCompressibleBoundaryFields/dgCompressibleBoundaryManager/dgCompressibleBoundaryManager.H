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
    Foam::dgCompressibleBoundaryManager

Description
    Manager for patch-level compressible DG boundary conditions.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCompressibleBoundaryManager_H
#define Foam_dgCompressibleBoundaryManager_H

#include "dictionary.H"
#include "IOdictionary.H"
#include "regIOobject.H"
#include "DynamicList.H"
#include "labelList.H"
#include "dgGeomMesh.H"
#include "dgCompressibleBoundaryField.H"
#include "dgBoundaryHelper.H"

namespace Foam
{

template<class Type> class dgField;
template<class Type> class GaussField;
class dgThermoConservative;

/*---------------------------------------------------------------------------*\
                Class dgCompressibleBoundaryManager Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgCompressibleBoundaryManager
 * \brief Manager for patch-wise compressible DG boundary conditions
 *
 * \details
 * Reads the compressible boundary-condition dictionary, constructs one
 * runtime-selected dgCompressibleBoundaryField for each supported physical
 * patch, and maintains the patch-to-condition lookup used during updates.
 *
 * The manager also stores the uniform internal primitive reference state and
 * uses it to initialize conservative DG fields before patch ghost states are
 * populated on boundary-face Gauss points.
 *
 * \note Compressible boundary conditions provide two weak forms:
 * weak Riemann ghost states, used by default for the governing equations in
 * DGFoam, and weak prescribed boundary values, reserved for auxiliary
 * equations solved with the LDG method.
 */
class dgCompressibleBoundaryManager
:
    public regIOobject
{
private:

    const dgGeomMesh& mesh_;
    const dgThermoConservative& thermo_;

    DynamicList<autoPtr<dgCompressibleBoundaryField>> bConditions_;
    labelList patchToCondition_;

    scalar internalP_;
    scalar internalT_;
    vector internalU_;

    void readInternalField(const dictionary& dict);
    void constructBoundaryConditions(const dictionary& dict);
    void resolveCouplings();

public:

    TypeName("dgCompressibleBoundaryManager");

    dgCompressibleBoundaryManager
    (
        const dictionary& boundaryDict,
        const dgThermoConservative& thermo,
        const dgGeomMesh& mesh
    );

    dgCompressibleBoundaryManager
    (
        const IOobject& io,
        const dgThermoConservative& thermo,
        const dgGeomMesh& mesh
    );

    dgCompressibleBoundaryManager
    (
        const dictionary& boundaryDict,
        const dgGeomMesh& mesh
    );

    dgCompressibleBoundaryManager
    (
        const IOobject& io,
        const dgGeomMesh& mesh
    );

    ~dgCompressibleBoundaryManager() = default;

    virtual bool writeData(Ostream& os) const override;

    inline label conditionIndex(const label patchID) const
    {
        return patchToCondition_[patchID];
    }

    label conditionIndex(const word& patchName) const;

    bool foundBoundaryField(const word& patchName) const;

    const dgCompressibleBoundaryField& boundaryField
    (
        const word& patchName
    ) const;

    dgCompressibleBoundaryField& boundaryFieldRef
    (
        const word& patchName
    );

    inline const scalar& internalP() const { return internalP_; }
    inline const scalar& internalT() const { return internalT_; }
    inline const vector& internalU() const { return internalU_; }

    void initializeConservatives
    (
        dgField<scalar>& rho,
        dgField<vector>& rhoU,
        dgField<scalar>& E
    ) const;

    void update
    (
        GaussField<scalar>& rhoG,
        GaussField<vector>& rhoUG,
        GaussField<scalar>& EG,
        GaussField<vector>& UG
    ) const;

    void correctSelfDiffusionFlux
    (
        GaussField<vector>& massDiffFlux
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
