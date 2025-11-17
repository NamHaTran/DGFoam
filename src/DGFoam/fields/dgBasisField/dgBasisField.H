/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::dgBasisField

Description
    Stores precomputed GaussField of all basis functions and their derivatives
    for a given cell in DG formulation.

SourceFiles
    dgBasisField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgBasisField_H
#define Foam_dgBasisField_H

#include "GaussField.H"
#include "dgGeomMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dgBasisField Declaration
\*---------------------------------------------------------------------------*/

class dgBasisField
{
private:

    //- Reference to geometric DG mesh
    const dgGeomMesh& mesh_;

    //- Cell ID associated with this basis field
    const label cellID_;

    //- Local face labels of the current cell
    labelList faceLabels_;

    //- Number of DoFs in this cell
    label nDof_;

    //- Basis GaussFields for all DoFs (scalar basis)
    List<GaussField<scalar>> basis_;

    //- Derivative basis GaussFields for all DoFs (vector basis)
    List<GaussField<vector>> dBasis_;

public:

    //- Construct from cell ID and geometric mesh
    dgBasisField(const label cellID, const dgGeomMesh& mesh);

    // Member Functions
        //- Compute scalar basis GaussFields (for all DoFs)
        void computeBasisField();

        //- Compute derivative basis GaussFields (for all DoFs)
        void computeDBasisField();

        // Access

        //- Return GaussField of the i-th basis function
        const GaussField<scalar>& getBasis(const label dof) const;

        //- Return GaussField of the i-th derivative basis function
        const GaussField<vector>& getDBasis(const label dof) const;

        //- Return number of DoF
        inline label nDof() const { return nDof_; }
};

} // End namespace Foam

// Include template-related inline implementation
#include "dgBasisFieldI.H"

#endif

// ************************************************************************* //
