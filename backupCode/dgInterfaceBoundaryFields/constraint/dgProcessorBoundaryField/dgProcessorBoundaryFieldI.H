/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DGFoam: Discontinuous Galerkin CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | GPU-friendly CFD solver framework
     \\/     M anipulation  |
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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dgProcessorBoundaryField

\*---------------------------------------------------------------------------*/

#ifndef dgProcessorBoundaryFieldI_H
#define dgProcessorBoundaryFieldI_H

#include "PstreamBuffers.H"
#include "Pstream.H"

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
    dgInterfaceBoundaryField<Type>(patch, dgMesh),
    procPatch_(refCast<const processorFvPatch>(patch)),
    neighbProc_(procPatch_.neighbProcNo())
{
    // Compute total Gauss points on this patch
    label totalGauss = 0;

    for (label i = 0; i < procPatch_.size(); ++i)
    {
        const label faceID = procPatch_.start() + i;
        const dgGeomFace* gFace = this->dgMesh_.faces()[faceID];

        totalGauss += gFace->connectivity().size();
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
    const int tag = Pstream::msgType();

    const label myProc = Pstream::myProcNo();

    // ----------------------------------
    // 1. Pack
    // ----------------------------------
    label k = 0;

    for (label i = 0; i < procPatch_.size(); ++i)
    {
        const label faceID = procPatch_.start() + i;
        const dgGeomFace* gFace = this->dgMesh_.faces()[faceID];

        const label ownerCellID =
            this->dgMesh_.getFvMesh().faceOwner()[faceID];

        const List<label>& conn = gFace->connectivity();
        const label nGauss = conn.size();

        for (label g = 0; g < nGauss; ++g)
        {
            sendBuf_[k] =
                field.gaussFields()[ownerCellID]
                     .faceField()
                     .minusValueOnGlobalFace(faceID, g);

            ++k;
        }
    }

    Pout<< "[Proc " << myProc << "] ";
    Pout << sendBuf_ << nl;

    // ----------------------------------
    // 2. MPI exchange
    // ----------------------------------

    {
        OPstream toNbr
        (
            Pstream::commsTypes::nonBlocking,
            tag,
            neighbProc_
        );

        toNbr << sendBuf_;
    }

    // Ensure sends are completed
    UPstream::waitRequests(-1);

    {
        IPstream fromNbr
        (
            Pstream::commsTypes::nonBlocking,
            tag,
            neighbProc_
        );

        fromNbr >> recvBuf_;
    }

    // ----------------------------------
    // 3. Unpack
    // ----------------------------------
    k = 0;

    for (label i = 0; i < procPatch_.size(); ++i)
    {
        const label faceID = procPatch_.start() + i;
        const dgGeomFace* gFace = this->dgMesh_.faces()[faceID];

        const label ownerCellID =
            this->dgMesh_.getFvMesh().faceOwner()[faceID];

        const List<label>& conn = gFace->connectivity();
        const label nGauss = conn.size();

        for (label g = 0; g < nGauss; ++g)
        {
            field.gaussFields()[ownerCellID]
                 .faceField()
                 .plusValueOnGlobalFace(faceID, conn[g]) = recvBuf_[k];

            ++k;
        }
    }
}

} // End namespace Foam

#endif
// ************************************************************************* //
