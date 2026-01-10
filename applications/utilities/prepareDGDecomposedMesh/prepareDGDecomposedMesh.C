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
    prepareDGDecomposedMesh

Group
    ultilities

Description
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "polyMesh.H"
#include "polyBoundaryMeshEntries.H"
#include "OSspecific.H"

// DG libs (must have)
#include "dgGeomMesh.H"
#include "dgField.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Polynomial order (temporary)
    const label pOrder = 1;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Creating DG Geometric Mesh with polynomial order "
        << pOrder << nl << endl;

    // ---------------------------------------------------------------------
    // Create DG geometric mesh (global mesh)
    // ---------------------------------------------------------------------

    dgGeomMesh dgMesh(mesh, pOrder);

    // ---------------------------------------------------------------------
    // Detect processor directories
    // ---------------------------------------------------------------------

    fileNameList caseDirs =
        Foam::readDir(runTime.path(), fileName::DIRECTORY);

    DynamicList<label> procIDsDyn;

    for (const fileName& dir : caseDirs)
    {
        if (dir.startsWith("processor"))
        {
            procIDsDyn.append(readLabel(dir.substr(9)));
        }
    }

    List<label> processorIDs(procIDsDyn);
    Foam::sort(processorIDs);

    if (processorIDs.empty())
    {
        FatalErrorInFunction
            << "No processor directories found." << nl
            << "Please run decomposePar before this tool."
            << exit(FatalError);
    }

    Info<< "Detected " << processorIDs.size()
        << " processor directories" << nl << endl;

    // ---------------------------------------------------------------------
    // Loop over processors
    // ---------------------------------------------------------------------

    forAll(processorIDs, pI)
    {
        const label procID = processorIDs[pI];

        Info<< "Processing processor " << procID << nl;

        const fileName polyMeshDir =
            runTime.path()
          / ("processor" + name(procID))
          / "constant/polyMesh";

        // --------------------------------------------------------------
        // Read addressing
        // --------------------------------------------------------------

        IOList<label> faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                polyMeshDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        IOList<label> cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                polyMeshDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read cell owner
        IOList<label> cellOwner
        (
            IOobject
            (
                "owner",
                polyMeshDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // --------------------------------------------------------------
        // Read boundary entries (dictionary-level, includes processor)
        // --------------------------------------------------------------

        IOobject boundaryIO
        (
            "boundary",
            polyMeshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        polyBoundaryMeshEntries boundaryEntries(boundaryIO);

        // --------------------------------------------------------------
        // Storage for DG face connectivity
        // Index: localProcFaceID
        // Value: Gauss connectivity permutation
        // --------------------------------------------------------------

        List<labelList> dgFaceConnectivity(faceProcAddressing.size());

        // --------------------------------------------------------------
        // Loop over boundary entries and process processor patches
        // --------------------------------------------------------------

        forAll(boundaryEntries, patchI)
        {
            const entry& e = boundaryEntries[patchI];
            const dictionary& dict = e.dict();

            const word patchType(dict.lookup("type"));

            if (patchType != "processor")
            {
                continue;
            }

            const label startFace =
                readLabel(dict.lookup("startFace"));

            const label nFaces =
                readLabel(dict.lookup("nFaces"));

            const label myProc =
                readLabel(dict.lookup("myProcNo"));

            //const label neighbProc =
            //    readLabel(dict.lookup("neighbProcNo"));

            /*
            Info<< "  Processor patch: " << e.keyword() << nl
                << "    myProc      = " << myProc << nl
                << "    neighbProc  = " << neighbProc << nl
                << "    startFace   = " << startFace << nl
                << "    nFaces      = " << nFaces << nl;
            */

            if (myProc != procID)
            {
                FatalErrorInFunction
                    << "Processor ID mismatch: boundary says myProc="
                    << myProc << " but directory is processor"
                    << procID
                    << exit(FatalError);
            }

            // ----------------------------------------------------------
            // Loop faces on this processor patch
            // ----------------------------------------------------------

            for (label i = 0; i < nFaces; ++i)
            {
                const label localProcFaceID = startFace + i;
                label globalFaceID =
                    faceProcAddressing[localProcFaceID];

                if (globalFaceID < 0)
                {
                    // Negative ID means face is owned by neighbor processor
                    // Must convert to positive ID before accessing global mesh
                    globalFaceID = -globalFaceID;
                }

                // Owner cell on processor mesh
                const label localProcCellID = cellOwner[localProcFaceID];

                const label globalCellID = cellProcAddressing[localProcCellID];

                // ------------------------------------------------------
                // Access DG face on global mesh
                // ------------------------------------------------------

                const dgGeomFace* gFace =
                    dgMesh.faces()[globalFaceID];

                const List<label>& connRef =
                    gFace->connectivity();

                List<label> procConnect(connRef.size(), -1);

                if (gFace->isOwner(globalCellID))
                {
                    // Owner (master) side
                    procConnect = connRef;
                }
                else
                {
                    // Neighbour (slave) side: reverse connectivity
                    forAll(connRef, k)
                    {
                        const label j = connRef[k];
                        procConnect[j] = k;
                    }
                }

                dgFaceConnectivity[localProcFaceID] = procConnect;
            }
        }

        // --------------------------------------------------------------
        // Write DG face connectivity
        // --------------------------------------------------------------

        fileName outFile =
            polyMeshDir / "dgFaceConnectivity";

        OFstream os(outFile);

        os  << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
            << "| =========                |                                                 |\n"
            << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
            << "|  \\    /   O peration     | Version:  2412                                  |\n"
            << "|   \\  /    A nd           | Website:  www.openfoam.com                      |\n"
            << "|    \\/     M anipulation  |                                                 |\n"
            << "\\*---------------------------------------------------------------------------*/\n"
            << "FoamFile\n"
            << "{\n"
            << "    version     2.0;\n"
            << "    format      ascii;\n"
            << "    class       labelListList;\n"
            << "    location    \"constant/polyMesh\";\n"
            << "    object      dgFaceConnectivity;\n"
            << "}\n"
            << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

        // Write data
        os << dgFaceConnectivity << nl;
    }

    Info<< nl << "prepareDGDecomposedMesh finished." << nl;

    return 0;
}


// ************************************************************************* //
