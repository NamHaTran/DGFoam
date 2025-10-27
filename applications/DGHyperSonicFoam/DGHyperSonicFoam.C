/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "dgRefFace.H"
#include "dgRefCell.H"
#include "dgGeomCell.H"
#include "dgGeomMesh.H"

// Test libs
#include "dgThermo.H"
#include "dgGeneralBoundaryManager.H"
//#include "cellDof.H"
#include "dofField.H"
//#include "cellGaussField.H"
//#include "faceGaussField.H"
#include "GaussField.H"
#include "dgGaussFieldLiteralScalarMath.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Declare pOrder
    const label pOrder = 2; // Polynomial order for basis functions

    #include "setRootCase.H"

    // *************************** OpenFOAM Initialization *************************** //
    // This includes the OpenFOAM initialization, which sets up the environment,
    // command line arguments, and the basic OpenFOAM data structures.
    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // Create fields
    //#include "createFields.H"

    // **************************** DGFoam Initialization **************************** //
    
    
    // Create the DG geometric mesh
    #include "createDGMesh.H"

    // Create the DG fields
    #include "createDGFields.H"

    // It's possible to iterate over every cell in a standard C++ for loop
    for (label cellI = 0; cellI < mesh.C().size(); cellI++)
    {
        // Create GaussField from DOF field
        GaussField<vector> UGauss
        (
            &UDof,
            cellI,
            dgMesh
        );

        if (cellI == 1000)
        {
            // Print results for cellID == 1000
            Info << "Cell " << cellI << ":\n";
            Info << "Original pGauss:\n";
            UGauss.print();
        }

        /*
        GaussField<scalar> pGauss
        (
            &pDof,
            cellI,
            dgMesh
        );

        // Interpolate values at Gauss points from DOF using basis functions
        UGauss.interpolateFromDof();
        pGauss.interpolateFromDof();

        GaussField<scalar> p1
        (
            cellI,
            dgMesh
        );
        GaussField<scalar> p2
        (
            cellI,
            dgMesh
        );
        GaussField<scalar> p3
        (
            cellI,
            dgMesh
        );
        GaussField<scalar> p4
        (
            cellI,
            dgMesh
        );
        GaussField<scalar> p5
        (
            cellI,
            dgMesh
        );
        GaussField<scalar> p6
        (
            cellI,
            dgMesh
        );

        // Test operators on GaussField
        p1 = pGauss;
        p2 = pGauss + 1e5;
        p3 = pGauss - 1e5;
        p4 = pGauss * 2.0;
        p5 = pGauss / 2.0;
        p6 = pow(pGauss, 2);

        if (cellI == 1000)
        {
            // Print results for cellID == 1000
            Info << "Cell " << cellI << ":\n";
            Info << "Original pGauss:\n";
            pGauss.print();
            Info << "p1 = pGauss:\n";
            p1.print();
            Info << "p2 = pGauss + 1e5:\n";
            p2.print();
            Info << "p3 = pGauss - 1e5:\n";
            p3.print();
            Info << "p4 = pGauss * 2.0:\n";
            p4.print();
            Info << "p5 = pGauss / 2.0:\n";
            p5.print();
            Info << "p6 = pow(pGauss, 2):\n";
            p6.print();
        }
        */
    }


    // Each cell is constructed of faces - these may either be internal or constitute a
    // boundary, or a patch in OpenFOAM terms; internal faces have an owner cell
    // and a neighbour.
    for (label faceI = 0; faceI < mesh.owner().size(); faceI++)
        if (faceI%40 == 0)
            Info << "Internal face " << faceI << " with centre at " << mesh.Cf()[faceI]
                 << " with owner cell " << mesh.owner()[faceI]
                 << " and neighbour " << mesh.neighbour()[faceI] << endl;
    Info << endl;

    // Boundary conditions may be accessed through the boundaryMesh object.
    // In reality, each boundary face is also included in the constant/polyMesh/faces
    // description. But, in that file, the internal faces are defined first.
    // In addition, the constant/polyMesh/boundary file defines the starting faceI
    // indices from which boundary face definitions start.
    // OpenFOAM also provides a macro definition for for loops over all entries
    // in a field or a list, which saves up on the amount of typing.
    forAll(mesh.boundaryMesh(), patchI)
        Info << "Patch " << patchI << ": " << mesh.boundary()[patchI].name() << " with "
             << mesh.boundary()[patchI].Cf().size() << " faces. Starts at total face "
             << mesh.boundary()[patchI].start() << endl;
    Info << endl;

    // Faces adjacent to boundaries may be accessed as follows.
    // Also, a useful thing to know about a face is its normal vector and face area.
    label patchFaceI(0);
    forAll(mesh.boundaryMesh(), patchI)
        Info << "Patch " << patchI << " has its face " << patchFaceI << " adjacent to cell "
             << mesh.boundary()[patchI].patch().faceCells()[patchFaceI]
             << ". It has normal vector " << mesh.boundary()[patchI].Sf()[patchFaceI]
             << " and surface area " << mag(mesh.boundary()[patchI].Sf()[patchFaceI])
             << endl;
    Info << endl;

    // For internal faces, method .Sf() can be called directly on the mesh instance.
    // Moreover, there is a shorthand method .magSf() which returns the surface area
    // as a scalar.
    // For internal faces, the normal vector points from the owner to the neighbour
    // and the owner has a smaller cellI index than the neighbour. For boundary faces,
    // the normals always point outside of the domain (they have "imaginary" neighbours
    // which do not exist).

    // It is possible to look at the points making up each face in more detail.
    // First, we define a few shorthands by getting references to the respective
    // objects in the mesh. These are defined as constants since we do not aim to
    // alter the mesh in any way.
    // NOTE: these lists refer to the physical definition of the mesh and thus
    // include boundary faces. Use can be made of the mesh.boundary()[patchI].Cf().size()
    // and mesh.boundary()[patchI].start() methods to check whether the face is internal
    // or lies on a boundary.
    const faceList& fcs = mesh.faces();
    const List<point>& pts = mesh.points();
    const List<point>& cents = mesh.faceCentres();

    forAll(fcs,faceI)
        if (faceI%80==0)
        {
            if (faceI<mesh.Cf().size())
                Info << "Internal face ";
            else
            {
                forAll(mesh.boundary(),patchI)
                    if ((mesh.boundary()[patchI].start()<= faceI) &&
                        (faceI < mesh.boundary()[patchI].start()+mesh.boundary()[patchI].Cf().size()))
                    {
                        Info << "Face on patch " << patchI << ", faceI ";
                        break; // exit the forAll loop prematurely
                    }
            }

            Info << faceI << " with centre at " << cents[faceI]
                 << " has " << fcs[faceI].size() << " vertices:";
            forAll(fcs[faceI],vertexI)
                // Note how fcs[faceI] holds the indices of points whose coordinates
                // are stored in the pts list.
                Info << " " << pts[fcs[faceI][vertexI]];
            Info << endl;
        }
    Info << endl;

    // In the original cavity tutorial, on which the test case is based,
    // the frontAndBack boundary is defined as and "empty" type. This is a special
    // BC case which may cause unexpected behaviour as its .Cf() field has size of 0.
    // Type of a patch may be checked to avoid running into this problem if there
    // is a substantial risk that an empty patch type will appear
    label patchID(0);
    const polyPatch& pp = mesh.boundaryMesh()[patchID];
    if (isA<emptyPolyPatch>(pp))
    {
        // patch patchID is of type "empty".
        Info << "You will not see this." << endl;
    }

    // Patches may also be retrieved from the mesh using their name. This could be
    // useful if the user were to refer to a particular patch from a dictionary
    // (like when you do when calculating forces on a particular patch).
    word patchName("movingWall");
    patchID = mesh.boundaryMesh().findPatchID(patchName);
    Info << "Retrieved patch " << patchName << " at index " << patchID << " using its name only." << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
