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
    //#include "createDGFields.H"

    // Read solver settings
    //#include "readSolverSettings.H"

    // TIME LOOP

        //#include "synchProcessors.H"

        // It's possible to iterate over every cell in a standard C++ for loop
        for (label cellI = 0; cellI < mesh.C().size(); cellI++)
        {
            dgMesh.cells()[cellI]->printDebugInfo();

            Foam::dgBasisField basisField(cellI, dgMesh);

            // -------------------------------- Update BCs --------------------------------- //
            // Declare primitive Gauss fields
            //GaussField<vector>& UG      = U.gaussFields()[cellI];
            //GaussField<scalar>& TG      = T.gaussFields()[cellI];
            //GaussField<scalar>& pG      = p.gaussFields()[cellI];
            //GaussField<scalar>& rhoG    = rho.gaussFields()[cellI];

            // Update ghost states
            //pBC     ->updateValue(pG);
            //TBC     ->updateValue(TG);
            //rhoBC   ->updateValue(rhoG);
            //UBC     ->updateValue(UG);

            // Update thermo
            //thermo->update(cellI);
            // ---------------------------------------------------------------------------- //


            // ---------------------------- Build equation terms -------------------------- //
            // Mass convection term
            /*
            GaussField<vector> rhoUG = rhoG*UG;
            dgGeneralPDETerm<scalar, vector> massConvTerm
            (
                "massConvTerm",
                dgFluxSolver::equationType::massTransport,
                cellI,
                dgMesh,
                rhoG,
                rhoUG,
                fluxSolverManager
            );

            massConvTerm.updateFlux();
            massConvTerm.calcResidual();

            // Get residual
            const List<scalar>& res = massConvTerm.R();

            Pout<< "Cell " << cellI << " residual: " << res << endl;
            */
        }

        //vectorInterfaceBC->update(U);

        // Reset necessary objects before next time step
        //#include "resetAll.H"
    
    // END TIME LOOP

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
