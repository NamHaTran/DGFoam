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

Description
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "dgGeomMesh.H"

// DG libs (must have)
#include "dgField.H"
#include "dgMath.H"
#include "dgGeneralBoundaryManager.H"
#include "dgProcessorBoundaryManager.H"
#include "dgFluxSolverManager.H"
#include "dgTimeDiscretization.H"
#include "dgThermoConservative.H"
#include "dgGeneralPDETerm.H"


// Test libs
#include "dgBasisField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Declare pOrder
    const label pOrder = 1; // Polynomial order for basis functions

    #include "setRootCase.H"

    // *************************** OpenFOAM Initialization *************************** //
    #include "createTime.H"
    #include "createMesh.H"

    // **************************** DGFoam Initialization **************************** //
    // Create the DG geometric mesh
    #include "createDGMesh.H"

    // Create the DG fields
    #include "createDGFields.H"

    // Read solver settings
    #include "readSolverSettings.H"

    // TIME LOOP
    while (runTime.run())
    {
        const scalar dtCandidate = timeDiscretization.deltaTValue();
        runTime.setDeltaT(dtCandidate, runTime.isAdjustTimeStep());

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl
            << "deltaT = " << runTime.deltaTValue() << nl << endl;

        while (timeDiscretization.loop())
        {
            #include "synchProcessors.H"
            #include "prepareStageFields.H"
            #include "assembleStageResiduals.H"
            #include "advanceStage.H"
        }

        timeDiscretization.reset();

        // Reset necessary objects before next time step
        #include "resetAll.H"
        runTime.write();
        runTime.printExecutionTime(Info);
    }
    
    // END TIME LOOP

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
