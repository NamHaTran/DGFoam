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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dgProcessorBoundaryField

Description
    Processor boundary handler for DG Gauss fields.

SourceFiles
    dgProcessorBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgProcessorBoundaryField_H
#define dgProcessorBoundaryField_H

#include "dgGeomMesh.H"
#include "processorFvPatch.H"
#include "dgField.H"
#include "Pstream.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dgProcessorBoundaryField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgProcessorBoundaryField
 * \brief DG processor boundary handler for parallel communication
 *
 * \tparam Type Field value type (e.g. scalar, vector, tensor)
 *
 * \details
 * Exchanges DG face Gauss-point plus-side values across processor patches.
 *
 * The communication uses explicit send/receive buffers and point-to-point
 * transfers on processorFvPatch.
 *
 * MPI tags are generated from processor ranks to avoid collisions.
 */
template<class Type>
class dgProcessorBoundaryField
{
private:

    /** Reference to processor fvPatch (type-safe) */
    const processorFvPatch& procPatch_;

    /** Reference to DG geometric mesh */
    const dgGeomMesh& dgMesh_;

    /** Neighbour processor rank */
    const label neighbProc_;

    /** Send buffer for Gauss-point values (mutable for const update) */
    mutable List<Type> sendBuf_;

    /** Receive buffer for Gauss-point values (mutable for const update) */
    mutable List<Type> recvBuf_;

public:

    // Constructors --------------------------------------------------------

    /**
     * \brief Construct from fvPatch and DG geometric mesh
     *
     * The input patch must be of type processor.
     *
     * \param patch  Finite-volume patch (processor)
     * \param dgMesh DG geometric mesh
     */
    dgProcessorBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh
    );

    //- Destructor
    virtual ~dgProcessorBoundaryField() = default;

    // Core communication --------------------------------------------------

    /**
     * \brief Exchange Gauss-point data across processor boundary
     *
     * Performs a send/receive cycle to update plus-side Gauss values on
     * processor faces.
     *
     * \param field DG field whose Gauss data will be exchanged
     */
    void update
    (
        dgField<Type>& field
    ) const;

    // Low-level communication helpers ------------------------------------

    /**
     * \brief Send all plus-side Gauss values to neighbour processor
     *
     * \param field DG field providing Gauss data
     */
    void sendAll
    (
        dgField<Type>& field
    ) const;

    /**
     * \brief Receive all plus-side Gauss values from neighbour processor
     *
     * \param field DG field to be updated with received data
     */
    void recvAll
    (
        dgField<Type>& field
    ) const;

    // MPI tag generation --------------------------------------------------

    /**
     * \brief Generate MPI send tag
     *
     * Tag scheme:
     * - tag = myProc*1e6 + neighbProc
     *
     * \return MPI tag for send
     */
    int genSendTag() const
    {
        const label myProc = Pstream::myProcNo();

        if (myProc >= 200000 || neighbProc_ >= 200000)
        {
            FatalErrorInFunction
                << "Processor rank too large for tag scheme: "
                << myProc << " " << neighbProc_
                << exit(FatalError);
        }

        return myProc*1e6 + neighbProc_;
    }

    /**
     * \brief Generate MPI receive tag
     *
     * Tag scheme:
     * - tag = neighbProc*1e6 + myProc
     *
     * \return MPI tag for receive
     */
    int genRecvTag() const
    {
        const label myProc = Pstream::myProcNo();

        if (myProc >= 200000 || neighbProc_ >= 200000)
        {
            FatalErrorInFunction
                << "Processor rank too large for tag scheme: "
                << myProc << " " << neighbProc_
                << exit(FatalError);
        }

        return neighbProc_*1e6 + myProc;
    }

    // Access functions ----------------------------------------------------

    /** Return processor patch */
    inline const processorFvPatch& patch() const { return procPatch_; }

    /** Return DG geometric mesh */
    inline const dgGeomMesh& dgMesh() const { return dgMesh_; }

    /** Return patch name */
    inline const word& name() const { return procPatch_.name(); }
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgProcessorBoundaryFieldI.H"

#endif
// ************************************************************************* //
