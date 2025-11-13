/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR by YOUR NAME
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

\*---------------------------------------------------------------------------*/

#include "faceGaussField.H"
#include "vector.H"

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Default constructor
template<class Type>
Foam::faceGaussField<Type>::faceGaussField()
:
    mesh_(nullptr),
    cellID_(-1),
    cell_(nullptr),
    facesID_(nullptr),
    nFaces_(0),
    nGaussPerFace_(0),
    nGauss_(0),
    faces_(),
    cellsDof_(),
    gaussOffset_(),
    plusValues_(),
    minusValues_()
{}

template<class Type>
Foam::faceGaussField<Type>::faceGaussField
(
    List<const cellDof<Type>*>& cellsDof,
    const dgGeomMesh* mesh
)
:
    mesh_(mesh),
    cellID_(cellsDof[0]->cellID()),
    cell_(mesh->cells()[cellID_]),
    facesID_(&(cell_->faces())),
    nFaces_(facesID_->size()),
    faces_(nFaces_),
    cellsDof_(cellsDof)
{
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        faces_[faceI] = mesh_->faces()[(*facesID_)[faceI]];
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    gaussOffset_.setSize(nFaces_);
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        gaussOffset_[faceI] = faceI * nGaussPerFace_;
    }

    plusValues_.setSize(nGauss_);
    minusValues_.setSize(nGauss_);

    interpolateFromDof();
}

template<class Type>
Foam::faceGaussField<Type>::faceGaussField
(
    const label cellID,
    const dgGeomMesh* mesh,
    const Type& initialValues
)
:
    mesh_(mesh),
    cellID_(cellID),
    cell_(mesh->cells()[cellID_]),
    facesID_(&(cell_->faces())),
    nFaces_(facesID_->size()),
    faces_(nFaces_),
    cellsDof_(0)
{
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        faces_[faceI] = mesh_->faces()[(*facesID_)[faceI]];
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    gaussOffset_.setSize(nFaces_);
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        gaussOffset_[faceI] = faceI * nGaussPerFace_;
    }

    plusValues_.setSize(nGauss_, initialValues);
    minusValues_.setSize(nGauss_, initialValues);
}


template<class Type>
Foam::faceGaussField<Type>::faceGaussField
(
    const label cellID,
    const dgGeomMesh* mesh
)
:
    mesh_(mesh),
    cellID_(cellID),
    cell_(mesh->cells()[cellID_]),
    facesID_(&(cell_->faces())),
    nFaces_(facesID_->size()),
    faces_(nFaces_),
    cellsDof_(0)
{
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        faces_[faceI] = mesh_->faces()[(*facesID_)[faceI]];
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    gaussOffset_.setSize(nFaces_);
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        gaussOffset_[faceI] = faceI * nGaussPerFace_;
    }

    plusValues_.setSize(nGauss_);
    minusValues_.setSize(nGauss_);
}


template<class Type>
Foam::faceGaussField<Type>::faceGaussField
(
    const faceGaussField<Type>& other
)
:
    // Copy mesh and geometry references
    mesh_(other.mesh_),             // Pointer to mesh (non-owning)
    cellID_(other.cellID_),         // ID of the associated cell
    cell_(other.cell_),             // Pointer to geometric cell (non-owning)
    facesID_(other.facesID_),       // Pointer to list of face IDs (non-owning)

    // Copy geometry-related information
    nFaces_(other.nFaces_),         // Number of faces in the cell
    nGaussPerFace_(other.nGaussPerFace_), // Number of Gauss points per face
    nGauss_(other.nGauss_),         // Total number of Gauss points on all faces

    // Copy pointer lists (shallow copy, non-owning)
    faces_(other.faces_),           // List of dgGeomFace* (non-owning)
    cellsDof_(other.cellsDof_),     // List of DOF pointers for owner and neighbour cells

    // Copy indexing and Gauss-point offset data
    gaussOffset_(other.gaussOffset_), // Per-face Gauss offset list

    // Deep copy of actual Gauss field values
    plusValues_(other.plusValues_),   // '+' side (neighbour) Gauss-point values
    minusValues_(other.minusValues_)  // '-' side (owner) Gauss-point values
{
    // No dynamic allocation here
    // All pointers are non-owning, safe for shallow copy
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
void Foam::faceGaussField<Type>::interpolateFromDof()
{
    const cellDof<Type>& cellPDof = *cellsDof_[0];

    for (label localFaceI = 0; localFaceI < nFaces_; ++localFaceI)
    {
        const label faceID = (*facesID_)[localFaceI];
        const label owner = mesh_->faceOwner()[faceID];
        const bool isOwner = (owner == cellID_);
        const bool isBoundaryFace = (faceID >= mesh_->nInternalFaces());

        const dgGeomFace& face = *faces_[localFaceI];

        // Get Gauss points for owner sides
        const List<vector>& gpOwner = face.gaussPointsOwner();

        // Get basis functions for owner sides
        const List<List<scalar>>& ownerBasis = face.ownerBasis();

        const label offset = gaussOffset_[localFaceI];

        if (isOwner)
        {   
            forAll(gpOwner, i)
            {
                minusValues_[offset + i] = pTraits<Type>::zero;
                for (label k = 0; k < cellPDof.nDof(); ++k)
                {
                    minusValues_[offset + i] += cellPDof[k] * ownerBasis[i][k];
                }
            }

            if (isBoundaryFace)
            {
                forAll(gpOwner, i)
                {
                    plusValues_[offset + i] = pTraits<Type>::max; //pTraits<Type>::zero;
                }
            }
            else
            {
                // Get Gauss points for neighbor sides
                const List<vector>& gpNeigh = face.gaussPointsNeighbor();

                // Get basis functions for neighbor sides
                const List<List<scalar>>& neighborBasis = face.neighborBasis();

                const cellDof<Type>& cellNDof = *cellsDof_[localFaceI + 1];
                forAll(gpNeigh, i)
                {
                    plusValues_[offset + i] = pTraits<Type>::zero;
                    for (label k = 0; k < cellNDof.nDof(); ++k)
                    {
                        plusValues_[offset + i] += cellNDof[k] * neighborBasis[i][k];
                    }
                }
            }
        }
        else
        {
            // This face has both owner and neighbor sides, so it is internal
            
            // Get Gauss points for neighbor sides
            const List<vector>& gpNeigh = face.gaussPointsNeighbor();

            // Get basis functions for neighbor sides
            const List<List<scalar>>& neighborBasis = face.neighborBasis();

            const cellDof<Type>& cellNDof = *cellsDof_[localFaceI + 1];
            
            forAll(gpNeigh, i)
            {
                minusValues_[offset + i] = pTraits<Type>::zero;
                for (label k = 0; k < cellPDof.nDof(); ++k)
                {
                    minusValues_[offset + i] += cellPDof[k] * neighborBasis[i][k];
                }
            }

            forAll(gpOwner, i)
            {
                plusValues_[offset + i] = pTraits<Type>::zero;
                for (label k = 0; k < cellNDof.nDof(); ++k)
                {
                    plusValues_[offset + i] += cellNDof[k] * ownerBasis[i][k];
                }
            }
        }
    }
}

template<class Type>
Foam::faceGaussField<Type>& Foam::faceGaussField<Type>::operator=
(
    const faceGaussField<Type>& other
)
{
    if (this == &other)
    {
        return *this;
    }

    // Copy geometric and topological references
    mesh_           = other.mesh_;          // Pointer to mesh (non-owning)
    cellID_         = other.cellID_;        // ID of the cell
    cell_           = other.cell_;          // Pointer to cell (non-owning)
    facesID_        = other.facesID_;       // Pointer to face IDs list (non-owning)

    // Copy geometry-related counts
    nFaces_         = other.nFaces_;        // Number of faces
    nGaussPerFace_  = other.nGaussPerFace_; // Number of Gauss points per face
    nGauss_         = other.nGauss_;        // Total number of Gauss points on all faces

    // Copy pointer lists (non-owning references)
    faces_          = other.faces_;         // List of dgGeomFace* (non-owning)
    cellsDof_       = other.cellsDof_;      // List of cellDof pointers (non-owning)

    // Copy indexing and offset information
    gaussOffset_    = other.gaussOffset_;   // Per-face offset of Gauss points

    // Deep copy of Gauss field values (these are actual data)
    plusValues_     = other.plusValues_;    // Values on the '+' (neighbor) side
    minusValues_    = other.minusValues_;   // Values on the '-' (owner) side

    return *this;
}

template<class Type>
Foam::faceGaussField<Type>& Foam::faceGaussField<Type>::operator=
(
    const Type& value
)
{
    for (label i = 0; i < minusValues_.size(); ++i)
    {
        minusValues_[i] = value;
        plusValues_[i] = value;
    }

    return *this;
}

// * * * * * * * * * * * * * * Template instantiations * * * * * * * * * * * * //

template class Foam::faceGaussField<Foam::scalar>;
template class Foam::faceGaussField<Foam::vector>;
template class Foam::faceGaussField<Foam::tensor>;
template class Foam::faceGaussField<Foam::symmTensor>;
template class Foam::faceGaussField<Foam::sphericalTensor>;

// ************************************************************************* //
} // End namespace Foam
