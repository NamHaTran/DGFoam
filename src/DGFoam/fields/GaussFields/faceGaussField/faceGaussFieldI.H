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

#ifndef Foam_faceGaussFieldI_H
#define Foam_faceGaussFieldI_H

#include "faceGaussField.H"

namespace Foam
{

// * * * * * * * * Constructors * * * * * * * * //

template<class Type>
faceGaussField<Type>::faceGaussField()
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


// From list of DOFs (P-cell + neighbour cells)
template<class Type>
faceGaussField<Type>::faceGaussField
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
    // Collect pointers to faces
    for (label f = 0; f < nFaces_; ++f)
    {
        faces_[f] = mesh_->faces()[(*facesID_)[f]];
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    gaussOffset_.setSize(nFaces_);
    for (label f = 0; f < nFaces_; ++f)
    {
        gaussOffset_[f] = f * nGaussPerFace_;
    }

    plusValues_.setSize(nGauss_);
    minusValues_.setSize(nGauss_);

    interpolateFromDof();
}


// Value-initialising constructor
template<class Type>
faceGaussField<Type>::faceGaussField
(
    const label cellID,
    const dgGeomMesh* mesh,
    const Type& initVal
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
    for (label f = 0; f < nFaces_; ++f)
    {
        faces_[f] = mesh_->faces()[(*facesID_)[f]];
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    gaussOffset_.setSize(nFaces_);
    for (label f = 0; f < nFaces_; ++f)
    {
        gaussOffset_[f] = f * nGaussPerFace_;
    }

    plusValues_.setSize(nGauss_, initVal);
    minusValues_.setSize(nGauss_, initVal);
}


// Zero-initialising constructor
template<class Type>
faceGaussField<Type>::faceGaussField
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
    for (label f = 0; f < nFaces_; ++f)
    {
        faces_[f] = mesh_->faces()[(*facesID_)[f]];
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    gaussOffset_.setSize(nFaces_);
    for (label f = 0; f < nFaces_; ++f)
    {
        gaussOffset_[f] = f * nGaussPerFace_;
    }

    plusValues_.setSize(nGauss_);
    minusValues_.setSize(nGauss_);
}


// Copy constructor
template<class Type>
faceGaussField<Type>::faceGaussField
(
    const faceGaussField<Type>& other
)
:
    mesh_(other.mesh_),
    cellID_(other.cellID_),
    cell_(other.cell_),
    facesID_(other.facesID_),
    nFaces_(other.nFaces_),
    nGaussPerFace_(other.nGaussPerFace_),
    nGauss_(other.nGauss_),
    faces_(other.faces_),
    cellsDof_(other.cellsDof_),
    gaussOffset_(other.gaussOffset_),
    plusValues_(other.plusValues_),
    minusValues_(other.minusValues_)
{}


// * * * * * * * * Member Functions * * * * * * * * //

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
                    plusValues_[offset + i] = minusValues_[offset + i]; //pTraits<Type>::max; //pTraits<Type>::zero;
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


// Assignment operator
template<class Type>
inline faceGaussField<Type>& faceGaussField<Type>::operator=
(
    const faceGaussField<Type>& other
)
{
    if (this == &other) return *this;

    mesh_          = other.mesh_;
    cellID_        = other.cellID_;
    cell_          = other.cell_;
    facesID_       = other.facesID_;
    nFaces_        = other.nFaces_;
    nGaussPerFace_ = other.nGaussPerFace_;
    nGauss_        = other.nGauss_;
    faces_         = other.faces_;
    cellsDof_      = other.cellsDof_;
    gaussOffset_   = other.gaussOffset_;
    plusValues_    = other.plusValues_;
    minusValues_   = other.minusValues_;

    return *this;
}


// Assign uniform value
template<class Type>
inline faceGaussField<Type>& faceGaussField<Type>::operator=
(
    const Type& v
)
{
    for (label i = 0; i < nGauss_; ++i)
    {
        plusValues_[i] = v; 
        minusValues_[i] = v;
    }
    return *this;
}

} // End namespace Foam

#endif

// ************************************************************************* //

