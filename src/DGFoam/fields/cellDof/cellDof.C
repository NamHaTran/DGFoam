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

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellDof_C
#define Foam_cellDof_C

#include "cellDof.H"

namespace Foam
{

template<class Type>
Foam::cellDof<Type>::cellDof
(
    const label cellID,
    dgGeomMesh& dgMesh,
    const UList<Type>& inputDof
)
:
    cellID_(cellID),
    dgMesh_(dgMesh),                    // Store reference to mesh
    cell_(*dgMesh.cells()[cellID]),     // Reference to geometric cell
    nDof_(0)
{
    // Get number of Dof from basis size
    nDof_ = cell_.basis().size();

    // Resize internal Dof list
    dof_.setSize(nDof_);

    const label inputSize = inputDof.size();
    const label minSize = min(inputSize, nDof_);

    // Copy values from inputDof to internal list
    for (label i = 0; i < minSize; ++i)
    {
        dof_[i] = inputDof[i];
    }

    // Fill remaining values with zero (default constructor)
    for (label i = minSize; i < nDof_; ++i)
    {
        dof_[i] = pTraits<Type>::zero;
    }
}

// * * * * * * * * * * * * * * Template instantiation * * * * * * * * * * * * //
template class Foam::cellDof<scalar>;
template class Foam::cellDof<vector>;
template class Foam::cellDof<tensor>;
template class Foam::cellDof<symmTensor>;
template class Foam::cellDof<sphericalTensor>;

} // End namespace Foam

#endif

// ************************************************************************* //
