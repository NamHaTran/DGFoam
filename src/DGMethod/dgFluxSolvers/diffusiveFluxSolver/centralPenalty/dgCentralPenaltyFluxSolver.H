/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::dgCentralPenaltyFluxSolver

Description
    Central diffusive flux with a CIP penalty term.

SourceFiles
    dgCentralPenaltyFluxSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCentralPenaltyFluxSolver_H
#define Foam_dgCentralPenaltyFluxSolver_H

#include "dgFluxSolver.H"
#include "dgField.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class dgCentralPenaltyFluxSolver Declaration
\*---------------------------------------------------------------------------*/

class dgCentralPenaltyFluxSolver
:
    public dgFluxSolver
{
private:

    /** Viscosity field used in the penalty coefficient */
    const dgField<scalar>& mu_;

    /** User coefficient in the CIP penalty */
    scalar CIP_;

    /** Name of the viscosity field in the object registry */
    word muName_;

    static word readMuName(const dictionary& dict);

    scalar penaltyCoeff
    (
        const label cellID,
        const label fI,
        const faceGaussField<scalar>& muFace
    ) const;

public:

    TypeName("centralPenalty");

    dgCentralPenaltyFluxSolver
    (
        const word& name,
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    virtual ~dgCentralPenaltyFluxSolver() = default;

    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<vector>& F,
        const faceGaussField<scalar>& U
    ) override;

    virtual void computeFlux
    (
        const label cellID,
        faceGaussField<tensor>& F,
        const faceGaussField<vector>& U
    ) override;

    virtual void reset() override {}

    virtual void read(const dictionary& dict) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
