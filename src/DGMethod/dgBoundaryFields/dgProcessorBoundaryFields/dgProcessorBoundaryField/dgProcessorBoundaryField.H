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

Description


SourceFiles
    dgProcessorBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgProcessorBoundaryField_H
#define dgProcessorBoundaryField_H

#include "dgGeomMesh.H"
#include "processorFvPatch.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dgProcessorBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgProcessorBoundaryField
{
private:

    //- Reference to processor patch (type-safe)
    const processorFvPatch& procPatch_;

    //- Reference to DG geometric mesh (for DG-level cell/face connectivity)
    const dgGeomMesh& dgMesh_;

    //- Neighbour processor rank
    const label neighbProc_;

    //- Communication buffers
    mutable List<Type> sendBuf_;
    mutable List<Type> recvBuf_;

public:

    //- Constructor from name and dictionary
    dgProcessorBoundaryField
    (
        const fvPatch& patch,
        const dgGeomMesh& dgMesh
    );

    //- Destructor
    virtual ~dgProcessorBoundaryField() = default;

    //- Exchange Gauss data across processor patch
    void update
    (
        dgField<Type>& field
    ) const;

    // Member Functions

        //- Send all plus values to neighbour processor
        void sendAll
        (
            dgField<Type>& field
        ) const;

        //- Receive all plus values from neighbour processor
        void recvAll
        (
            dgField<Type>& field
        ) const;

        //- Generate tag for MPI communication
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

        //- Check patch type
        /*
        void checkPatchType() const
        {
            if (
                this->procPatch_.type() != "processor"
            )
            {
                FatalErrorInFunction
                    << "Boundary condition " << type() <<" can only "
                    << "be applied to patch type: \n"
                    << "    " << "processor \n"
                    << "but patch "
                    << this->procPatch_.name() << " is of type " << this->procPatch_.type() << nl
                    << exit(FatalError);
            }
        }*/

    // Access Functions

        inline const processorFvPatch& patch() const { return procPatch_; }

        inline const dgGeomMesh& dgMesh() const { return dgMesh_; }

        inline const word& name() const { return procPatch_.name(); }
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgProcessorBoundaryFieldI.H"

#endif
// ************************************************************************* //
