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
    Foam::dgProcessorBoundaryField

\*---------------------------------------------------------------------------*/

#ifndef dgProcessorBoundaryFieldI_H
#define dgProcessorBoundaryFieldI_H

#include "PstreamBuffers.H"
#include "Pstream.H"
#include "UOPstream.H"
#include "UIPstream.H"
#include "transform.H"

namespace Foam
{

// * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

template<class Type>
dgProcessorBoundaryField<Type>::dgProcessorBoundaryField
(
    const fvPatch& patch,
    const dgGeomMesh& dgMesh
)
:
    procPatch_(refCast<const processorFvPatch>(patch)),
    dgMesh_(dgMesh),
    neighbProc_(procPatch_.neighbProcNo())
{
    // Compute total Gauss points on this patch
    label totalGauss = 0;

    for (label i = 0; i < procPatch_.size(); ++i)
    {
        const label faceID = procPatch_.start() + i;
        const dgGeomFace* gFace = this->dgMesh_.faces()[faceID];

        totalGauss += gFace->nGauss();
    }

    // Allocate send and receive buffers
    sendBuf_.setSize(totalGauss);
    recvBuf_.setSize(totalGauss);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class Type>
void dgProcessorBoundaryField<Type>::update
(
    dgField<Type>& field
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    // Reserve one temporary message tag for this entire processor exchange.
    // This isolates the DG communication from any other MPI traffic that may
    // happen in the solver during the same stage.
    const int oldTag = UPstream::incrMsgType();

    // PstreamBuffers with nonBlocking comms implements the common OpenFOAM
    // two-phase pattern:
    // 1. stream all outgoing data into per-neighbour buffers,
    // 2. call finishedSends() to exchange the buffered messages,
    // 3. read the incoming buffers back in the same serialization order.
    //
    // This avoids the deadlock-prone "both ranks do blocking send first,
    // then blocking receive" pattern.
    PstreamBuffers pBufs
    (
        UPstream::commsTypes::nonBlocking,
        UPstream::msgType()
    );

    sendAll(field, pBufs);
    pBufs.finishedSends();
    recvAll(field, pBufs);

    UPstream::msgType(oldTag);
}


template<class Type>
void dgProcessorBoundaryField<Type>::sendAll
(
    dgField<Type>& field,
    PstreamBuffers& pBufs
) const
{
    // If the field owns modal DoFs, the message has two serialized sections:
    // 1. all face Gauss-point minus states on this processor patch,
    // 2. one packed per-face DoF payload.
    //
    // recvAll() must deserialize in exactly the same order.
    const bool hasDof = field.hasDof();

    // ----------------------------------
    // Pack
    // ----------------------------------
    label k = 0;

    for (label i = 0; i < procPatch_.size(); ++i)
    {
        const label faceID = procPatch_.start() + i;
        const dgGeomFace* gFace = this->dgMesh_.faces()[faceID];

        const label ownerCellID =
            this->dgMesh_.getFvMesh().faceOwner()[faceID];

        const label nGauss = gFace->nGauss();

        for (label g = 0; g < nGauss; ++g)
        {
            // For a processor face we export the local owner-side trace
            // ("minus" state from this rank's point of view). The neighbour
            // will reconstruct its local plus-side data from this payload.
            sendBuf_[k] =
                field.gaussFields()[ownerCellID]
                     .faceField()
                     .minusValueOnGlobalFace(faceID, g);

            ++k;
        }
    }

    if (hasDof)
    {
        labelList faceDofSizes(procPatch_.size(), Zero);

        for (label i = 0; i < procPatch_.size(); ++i)
        {
            const label faceID = procPatch_.start() + i;
            const label ownerCellID =
                this->dgMesh_.getFvMesh().faceOwner()[faceID];

            faceDofSizes[i] = field.dof()[ownerCellID].nDof();
        }

        sendDofBuf_.resize(faceDofSizes);

        for (label i = 0; i < procPatch_.size(); ++i)
        {
            const label faceID = procPatch_.start() + i;
            const label ownerCellID =
                this->dgMesh_.getFvMesh().faceOwner()[faceID];
            const cellDof<Type>& ownerDof = field.dof()[ownerCellID];
            SubList<Type> faceDofs = sendDofBuf_[i];

            for (label dofI = 0; dofI < ownerDof.nDof(); ++dofI)
            {
                faceDofs[dofI] = ownerDof[dofI];
            }
        }
    }

    // ----------------------------------
    // Send to neighbour
    // ----------------------------------

    {
        // UOPstream(neighbProc_, pBufs) does not send immediately with a
        // blocking MPI_Send. It serializes into the neighbour-specific buffer
        // held by pBufs. The actual exchange is triggered later by
        // pBufs.finishedSends() in update()/manager::update().
        UOPstream toNbr(neighbProc_, pBufs);

        toNbr << sendBuf_;

        if (hasDof)
        {
            toNbr << sendDofBuf_;
        }
    }
}

template<class Type>
void dgProcessorBoundaryField<Type>::recvAll
(
    dgField<Type>& field,
    PstreamBuffers& pBufs
) const
{
    // recvAll() must consume the stream in exactly the order produced by
    // sendAll():
    // 1. the packed Gauss-point values list,
    // 2. optionally one packed DoF payload for all faces.
    const bool hasDof = field.hasDof();

    // ----------------------------------
    // Receive from neighbour
    // ----------------------------------

    // UIPstream reads from the neighbour-specific receive buffer that was
    // populated by pBufs.finishedSends(). No extra point-to-point MPI call
    // happens here.
    UIPstream fromNbr(neighbProc_, pBufs);

    fromNbr >> recvBuf_;

    if (hasDof)
    {
        fromNbr >> recvDofBuf_;
    }

    // ----------------------------------
    // Unpack
    // ----------------------------------
    label k = 0;

    for (label i = 0; i < procPatch_.size(); ++i)
    {
        const label faceID = procPatch_.start() + i;
        const dgGeomFace* gFace = this->dgMesh_.faces()[faceID];

        const label ownerCellID =
            this->dgMesh_.getFvMesh().faceOwner()[faceID];

        const label nGauss = gFace->nGauss();
        faceGaussField<Type>& faceField =
            field.gaussFields()[ownerCellID].faceField();

        for (label g = 0; g < nGauss; ++g)
        {
            // The neighbour's exported minus state becomes this rank's local
            // plus state on the same global face.
            faceField.plusValueOnGlobalFace(faceID, g)
                = recvBuf_[k];

            ++k;
        }

        if (hasDof)
        {
            List<Type> recvDof
            (
                recvDofBuf_.localSize(i),
                pTraits<Type>::zero
            );
            const SubList<Type> faceRecvDof = recvDofBuf_[i];

            for (label dofI = 0; dofI < recvDof.size(); ++dofI)
            {
                recvDof[dofI] = faceRecvDof[dofI];
            }

            // Store a copy in face-local buffered storage so later limiter /
            // reconstruction code can query neighbourCellDof() safely without
            // depending on the communication stream lifetime.
            faceField.setBufferedNeighbourDof
            (
                faceField.localFaceIndex(faceID),
                cellDof<Type>(-1, recvDof.size(), recvDof)
            );
        }
    }
}

} // End namespace Foam

#endif
// ************************************************************************* //
