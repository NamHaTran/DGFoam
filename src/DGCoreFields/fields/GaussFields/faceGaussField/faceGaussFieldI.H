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
    Foam::faceGaussField<Type>

Description
    Inline implementations for Foam::faceGaussField.

SourceFiles
    faceGaussFieldI.H

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
    bufferedNbrDof_(),
    bufferedNbrDofValid_(),
    gaussOffset_(),
    ownerPatchIDs_(),
    n_(),
    faceAreas_(),
    plusValues_(),
    minusValues_(),
    midValues_(),
    flux_()
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
    cellsDof_(cellsDof),
    bufferedNbrDof_(nFaces_),
    bufferedNbrDofValid_(nFaces_, false),
    n_(nFaces_),
    faceAreas_(nFaces_)
{
    // Collect pointers to faces
    for (label f = 0; f < nFaces_; ++f)
    {
        const label faceID = (*facesID_)[f];
        faces_[f] = mesh_->faces()[faceID];

        const vector faceNormal = faces_[f]->normal();
        n_[f] =
            (mesh_->faceOwner()[faceID] == cellID_ ? faceNormal : -faceNormal);

        faceAreas_[f] = faces_[f]->area();
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    // Get Gauss offset and face's owner patch's ID
    gaussOffset_.setSize(nFaces_);
    ownerPatchIDs_.setSize(nFaces_);
    for (label f = 0; f < nFaces_; ++f)
    {
        gaussOffset_[f] = f * nGaussPerFace_;

        if (faces_[f]->isBoundary())
        {
            ownerPatchIDs_[f] = faces_[f]->patchID();
        }
        else
        {
            ownerPatchIDs_[f] = -1;
        }
    }


    plusValues_.setSize(nGauss_);
    minusValues_.setSize(nGauss_);
    midValues_.setSize(nGauss_);

    interpolateFromDof();
    midValues_ = plusValues_;
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
    cellsDof_(0),
    bufferedNbrDof_(nFaces_),
    bufferedNbrDofValid_(nFaces_, false),
    n_(nFaces_),
    faceAreas_(nFaces_)
{
    for (label f = 0; f < nFaces_; ++f)
    {
        const label faceID = (*facesID_)[f];
        faces_[f] = mesh_->faces()[faceID];

        const vector faceNormal = faces_[f]->normal();
        n_[f] =
            (mesh_->faceOwner()[faceID] == cellID_ ? faceNormal : -faceNormal);
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    // Get Gauss offset and face's owner patch's ID
    gaussOffset_.setSize(nFaces_);
    ownerPatchIDs_.setSize(nFaces_);
    for (label f = 0; f < nFaces_; ++f)
    {
        gaussOffset_[f] = f * nGaussPerFace_;

        if (faces_[f]->isBoundary())
        {
            ownerPatchIDs_[f] = faces_[f]->patchID();
        }
        else
        {
            ownerPatchIDs_[f] = -1;
        }
    }

    plusValues_.setSize(nGauss_, initVal);
    minusValues_.setSize(nGauss_, initVal);
    midValues_.setSize(nGauss_, initVal);
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
    cellsDof_(0),
    bufferedNbrDof_(nFaces_),
    bufferedNbrDofValid_(nFaces_, false),
    n_(nFaces_),
    faceAreas_(nFaces_)
{
    for (label f = 0; f < nFaces_; ++f)
    {
        const label faceID = (*facesID_)[f];
        faces_[f] = mesh_->faces()[faceID];

        const vector faceNormal = faces_[f]->normal();
        n_[f] =
            (mesh_->faceOwner()[faceID] == cellID_ ? faceNormal : -faceNormal);
    }

    nGaussPerFace_ = faces_[0]->gaussPointsOwner().size();
    nGauss_ = nFaces_ * nGaussPerFace_;

    // Get Gauss offset and face's owner patch's ID
    gaussOffset_.setSize(nFaces_);
    ownerPatchIDs_.setSize(nFaces_);
    for (label f = 0; f < nFaces_; ++f)
    {
        gaussOffset_[f] = f * nGaussPerFace_;

        if (faces_[f]->isBoundary())
        {
            ownerPatchIDs_[f] = faces_[f]->patchID();
        }
        else
        {
            ownerPatchIDs_[f] = -1;
        }
    }

    plusValues_.setSize(nGauss_);
    minusValues_.setSize(nGauss_);
    midValues_.setSize(nGauss_, pTraits<Type>::zero);
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
    bufferedNbrDof_(other.bufferedNbrDof_),
    bufferedNbrDofValid_(other.bufferedNbrDofValid_),
    gaussOffset_(other.gaussOffset_),
    ownerPatchIDs_(other.ownerPatchIDs_),
    n_(other.n_),
    faceAreas_(other.faceAreas_),
    plusValues_(other.plusValues_),
    minusValues_(other.minusValues_),
    midValues_(other.midValues_),
    flux_(other.flux_)
{}

// Move constructor
template<class Type>
faceGaussField<Type>::faceGaussField
(
    faceGaussField<Type>&& other
) noexcept
:
    mesh_(other.mesh_),
    cellID_(other.cellID_),
    cell_(other.cell_),
    facesID_(other.facesID_),
    nFaces_(other.nFaces_),
    nGaussPerFace_(other.nGaussPerFace_),
    nGauss_(other.nGauss_),
    faces_(std::move(other.faces_)),
    cellsDof_(std::move(other.cellsDof_)),
    bufferedNbrDof_(std::move(other.bufferedNbrDof_)),
    bufferedNbrDofValid_(std::move(other.bufferedNbrDofValid_)),
    gaussOffset_(std::move(other.gaussOffset_)),
    ownerPatchIDs_(std::move(other.ownerPatchIDs_)),
    n_(std::move(other.n_)),
    faceAreas_(std::move(other.faceAreas_)),
    plusValues_(std::move(other.plusValues_)),
    minusValues_(std::move(other.minusValues_)),
    midValues_(std::move(other.midValues_)),
    flux_(std::move(other.flux_))
{}

// Construct from tmp
/*
template<class Type>
faceGaussField<Type>::faceGaussField
(
    const tmp<faceGaussField<Type>>& t
)
{
    if (t.movable())
    {
        this->operator=(std::move(const_cast<faceGaussField<Type>&>(t())));
        t.clear();
    }
    else
    {
        *this = t();
    }
}


template<class Type>
Foam::faceGaussField<Type>::faceGaussField
(
    const tmp<faceGaussField<Type>>& t
)
:
    faceGaussField<Type>(t())
{
    t.clear();
}
*/

template<class Type>
Foam::faceGaussField<Type>::faceGaussField
(
    const tmp<faceGaussField<Type>>& t
)
{
    if (t.movable())
    {
        this->operator=(std::move(const_cast<faceGaussField<Type>&>(t())));
        t.clear();
    }
    else
    {
        *this = t();
    }
}

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
        const bool hasNeighbourDof =
        (
            cellsDof_.size() > localFaceI + 1
         && cellsDof_[localFaceI + 1] != nullptr
        );
        bool isBoundaryFace = faces_[localFaceI]->isBoundary() || !hasNeighbourDof;

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

                const cellDof<Type>& cellNDof = neighbourCellDof(localFaceI);
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
            if (!hasNeighbourDof)
            {
                forAll(gpOwner, i)
                {
                    minusValues_[offset + i] = pTraits<Type>::zero;
                    for (label k = 0; k < cellPDof.nDof(); ++k)
                    {
                        minusValues_[offset + i] += cellPDof[k]*ownerBasis[i][k];
                    }

                    plusValues_[offset + i] = minusValues_[offset + i];
                }

                continue;
            }
            
            // Get Gauss points for neighbor sides
            const List<vector>& gpNeigh = face.gaussPointsNeighbor();

            // Get basis functions for neighbor sides
            const List<List<scalar>>& neighborBasis = face.neighborBasis();

            const cellDof<Type>& cellNDof = neighbourCellDof(localFaceI);
            
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
boundaryGaussField<Type> Foam::faceGaussField<Type>::extractBCGhostState() const
{
    label nBoundaryGauss = 0;

    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (isBoundary(faceI))
        {
            nBoundaryGauss += nGaussPerFace_;
        }
    }

    boundaryGaussField<Type> bgf(nBoundaryGauss);

    label bgi = 0;

    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (!isBoundary(faceI))
        {
            continue;
        }

        for (label gI = 0; gI < nGaussPerFace_; ++gI)
        {
            bgf[bgi++] = plusValueOnFace(faceI, gI);
        }
    }

    return bgf;
}


template<class Type>
void Foam::faceGaussField<Type>::assignBCGhostState
(
    const boundaryGaussField<Type>& bgf
)
{
    label expectedSize = 0;

    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (isBoundary(faceI))
        {
            expectedSize += nGaussPerFace_;
        }
    }

    if (bgf.size() != expectedSize)
    {
        FatalErrorInFunction
            << "boundaryGaussField size mismatch for cell " << cellID_ << nl
            << "Expected " << expectedSize << " values but received "
            << bgf.size() << abort(FatalError);
    }

    label bgi = 0;

    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (!isBoundary(faceI))
        {
            continue;
        }

        for (label gI = 0; gI < nGaussPerFace_; ++gI)
        {
            plusValueOnFace(faceI, gI) = bgf[bgi++];
        }
    }
}


template<class Type>
boundaryGaussField<Type> Foam::faceGaussField<Type>::extractMidValue() const
{
    label nBoundaryGauss = 0;

    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (isBoundary(faceI))
        {
            nBoundaryGauss += nGaussPerFace_;
        }
    }

    boundaryGaussField<Type> bgf(nBoundaryGauss);

    label bgi = 0;

    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (!isBoundary(faceI))
        {
            continue;
        }

        for (label gI = 0; gI < nGaussPerFace_; ++gI)
        {
            bgf[bgi++] = midValueOnFace(faceI, gI);
        }
    }

    return bgf;
}


template<class Type>
void Foam::faceGaussField<Type>::assignMidValue
(
    const boundaryGaussField<Type>& bgf
)
{
    label expectedSize = 0;

    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (isBoundary(faceI))
        {
            expectedSize += nGaussPerFace_;
        }
    }

    if (bgf.size() != expectedSize)
    {
        FatalErrorInFunction
            << "boundaryGaussField size mismatch for cell " << cellID_ << nl
            << "Expected " << expectedSize << " values but received "
            << bgf.size() << abort(FatalError);
    }

    label bgi = 0;

    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        if (!isBoundary(faceI))
        {
            continue;
        }

        for (label gI = 0; gI < nGaussPerFace_; ++gI)
        {
            midValueOnFace(faceI, gI) = bgf[bgi++];
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
    bufferedNbrDof_ = other.bufferedNbrDof_;
    bufferedNbrDofValid_ = other.bufferedNbrDofValid_;
    gaussOffset_   = other.gaussOffset_;
    ownerPatchIDs_ = other.ownerPatchIDs_;
    n_             = other.n_;
    faceAreas_     = other.faceAreas_;
    plusValues_    = other.plusValues_;
    minusValues_   = other.minusValues_;
    midValues_     = other.midValues_;
    flux_          = other.flux_;

    return *this;
}

// Move assignment
template<class Type>
inline faceGaussField<Type>& faceGaussField<Type>::operator=
(
    faceGaussField<Type>&& other
) noexcept
{
    if (this != &other)
    {
        mesh_          = other.mesh_;
        cellID_        = other.cellID_;
        cell_          = other.cell_;
        facesID_       = other.facesID_;
        nFaces_        = other.nFaces_;
        nGaussPerFace_ = other.nGaussPerFace_;
        nGauss_        = other.nGauss_;

        faces_         = std::move(other.faces_);
        cellsDof_      = std::move(other.cellsDof_);
        bufferedNbrDof_ = std::move(other.bufferedNbrDof_);
        bufferedNbrDofValid_ = std::move(other.bufferedNbrDofValid_);
        gaussOffset_   = std::move(other.gaussOffset_);
        ownerPatchIDs_ = std::move(other.ownerPatchIDs_);
        n_             = std::move(other.n_);
        faceAreas_     = std::move(other.faceAreas_);
        plusValues_    = std::move(other.plusValues_);
        minusValues_   = std::move(other.minusValues_);
        midValues_     = std::move(other.midValues_);
        flux_          = std::move(other.flux_);
    }

    return *this;
}

// Assign tmp
template<class Type>
inline faceGaussField<Type>& faceGaussField<Type>::operator=
(
    const tmp<faceGaussField<Type>>& t
)
{
    if (t.movable())
    {
        this->operator=(std::move(const_cast<faceGaussField<Type>&>(t())));
        t.clear();
    }
    else
    {
        *this = t();
    }

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
        midValues_[i] = v;
    }
    return *this;
}

} // End namespace Foam

#endif

// ************************************************************************* //
