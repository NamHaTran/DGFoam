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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::faceGaussField<Type>::faceGaussField
(
    const List<cellDof<Type>*>& cellsDof,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh),
    cellID_(cellsDof[0]->cellID()),
    cell_(*mesh.cells()[cellID_]),
    facesID_(cell_.faces()),
    nFaces_(facesID_.size()),
    faces_(nFaces_),
    cellsDof_(cellsDof),
    ctxPtr_(nullptr)
{
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        faces_[faceI] = mesh_.faces()[facesID_[faceI]];
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
    const List<cellDof<Type>*>& cellsDof,
    const GeometricField<Type, fvPatchField, volMesh>& foamField,
    const dgGeomMesh& mesh
)
:
    mesh_(mesh),
    cellID_(cellsDof[0]->cellID()),
    cell_(*mesh.cells()[cellID_]),
    facesID_(cell_.faces()),
    nFaces_(facesID_.size()),
    faces_(nFaces_),
    cellsDof_(cellsDof),
    ctxPtr_(nullptr)
{
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        faces_[faceI] = mesh_.faces()[facesID_[faceI]];
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
    const List<cellDof<Type>*>& cellsDof,
    const dgGeomMesh& mesh,
    const Type& initialValues
)
:
    mesh_(mesh),
    cellID_(cellsDof[0]->cellID()),
    cell_(*mesh.cells()[cellID_]),
    facesID_(cell_.faces()),
    nFaces_(facesID_.size()),
    faces_(nFaces_),
    cellsDof_(cellsDof),
    ctxPtr_(nullptr)
{
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        faces_[faceI] = mesh_.faces()[facesID_[faceI]];
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
    const dgGeomMesh& mesh,
    const Type& initialValues
)
:
    mesh_(mesh),
    cellID_(cellID),
    cell_(*mesh.cells()[cellID]),
    facesID_(cell_.faces()),
    nFaces_(facesID_.size()),
    faces_(nFaces_),
    cellsDof_(0),
    ctxPtr_(nullptr)
{
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        faces_[faceI] = mesh_.faces()[facesID_[faceI]];
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
    const dgGeomMesh& mesh
)
:
    mesh_(mesh),
    cellID_(cellID),
    cell_(*mesh.cells()[cellID]),
    facesID_(cell_.faces()),
    nFaces_(facesID_.size()),
    faces_(nFaces_),
    cellsDof_(0),
    ctxPtr_(nullptr)
{
    for (label faceI = 0; faceI < nFaces_; ++faceI)
    {
        faces_[faceI] = mesh_.faces()[facesID_[faceI]];
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
    minusValues_(other.minusValues_),
    ctxPtr_(other.ctxPtr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
void Foam::faceGaussField<Type>::interpolateFromDof()
{
    const cellDof<Type>& cellPDof = *cellsDof_[0];

    for (label localFaceI = 0; localFaceI < nFaces_; ++localFaceI)
    {
        const label faceID = facesID_[localFaceI];
        const label owner = mesh_.faceOwner()[faceID];
        const bool isOwner = (owner == cellID_);
        const bool isBoundaryFace = (faceID >= mesh_.nInternalFaces());

        const dgGeomFace& face = *faces_[localFaceI];

        // Get Gauss points for neighbor sides
        const List<vector>& gpOwner = face.gaussPointsOwner();
        const List<vector>& gpNeigh = face.gaussPointsNeighbor();

        // Get basis functions for neighbor sides
        const List<List<scalar>>& ownerBasis = face.ownerBasis();
        const List<List<scalar>>& neighborBasis = face.neighborBasis();

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
                forAll(gpNeigh, i)
                {
                    plusValues_[offset + i] = pTraits<Type>::zero;
                }
            }
        }
        else
        {
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
}

// * * * * * * * * * * * * * * Template instantiations * * * * * * * * * * * * //

template class Foam::faceGaussField<Foam::scalar>;
template class Foam::faceGaussField<Foam::vector>;
template class Foam::faceGaussField<Foam::tensor>;
template class Foam::faceGaussField<Foam::symmTensor>;
template class Foam::faceGaussField<Foam::sphericalTensor>;

// ************************************************************************* //

