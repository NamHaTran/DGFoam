/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Class
    Foam::dgGeomMesh

Description
    Geometric mesh class for DG.

    - Store all geometric faces and cells
    - Handle DG-related geometric construction from fvMesh
    - Share reference face/cell info via dgRefFace and dgRefCell

See also

SourceFiles
    dgGeomMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeomMesh_H
#define Foam_dgGeomMesh_H

#include "fvMesh.H"
#include "dgRefFace.H"
#include "dgGeomFace.H"
#include "dgGeomCell.H"
#include "dgRefCell.H"
#include "dgFacePosition.H"
#include "dgFaceType.H"
#include "dgCellType.H"

#include <memory>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dgGeomMesh Declaration
\*---------------------------------------------------------------------------*/

class dgGeomMesh
{
private:

    //- Reference to fvMesh
    const fvMesh& mesh_;

    //- Polynomial order
    label pOrder_;

    //- Reference face (shared among all faces)
    std::shared_ptr<dgRefFace> refFace_;

    //- Reference cells for different types
    std::shared_ptr<dgRefCell> refCellTet_;
    std::shared_ptr<dgRefCell> refCellHex_;
    std::shared_ptr<dgRefCell> refCellPrism_;
    std::shared_ptr<dgRefCell> refCellPyramid_;

    //- List of DG faces
    List<dgGeomFace*> faces_;

    //- List of DG cells
    List<dgGeomCell*> cells_;

    //- List of boundary faces
    labelList boundaryFaces_;

public:

    // Constructors

        //- Construct from fvMesh and polynomial order
        dgGeomMesh(const fvMesh& mesh, const label pOrder);


    //- Destructor
    ~dgGeomMesh();


    // Member Functions

        //- Return access to face list
        inline const List<dgGeomFace*>& faces() const
        {
            return faces_;
        }

        //- Return access to cell list
        inline const List<dgGeomCell*>& cells() const
        {
            return cells_;
        }

        //- Return polynomial order
        inline label pOrder() const
        {
            return pOrder_;
        }

        //- Return number of faces
        inline label nFaces() const
        {
            return faces_.size();
        }

        //- Return number of cells
        inline label nCells() const
        {
            return cells_.size();
        }

        //- Return face owner list from fvMesh
        inline const labelList& faceOwner() const { return mesh_.faceOwner(); }

        //- Return number of internal faces from fvMesh
        inline label nInternalFaces() const { return mesh_.nInternalFaces(); }

        //- Get boundary faces
        void getBoundaryFaces();

        //- Return boundary faces
        const Foam::labelList& boundaryFaces() const;

        //- Get patch id
        Foam::label getPatchID(const label faceID) const;

        //- Get local face ID within a patch
        Foam::label getLocalFaceID
        (
            const label faceID,
            const label patchID
        ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
