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
    Minimal density-based DG solver for inviscid/compressible experiments.

    Tutorial execution flow:
    1. Build the OpenFOAM time database and fvMesh.
    2. Wrap the fvMesh into a dgGeomMesh that stores DG geometry/basis data.
    3. Create conservative and primitive DG fields plus the thermo model.
    4. Create flux-solver and time-discretization managers from dgSchemes.
    5. For each physical time step:
       - compute a global explicit deltaT from the current DG state,
       - iterate over RK/Euler stages managed by dgTimeDiscretization,
       - prepare stage input states and boundary data,
       - assemble cell-local DG residuals,
       - advance the registered conservative fields.
    6. Write the updated state and proceed to the next physical step.

    Core dgFoam ideas used by this tutorial:
    - dgGeomMesh extends fvMesh with DG geometry, quadrature, basis, and
      cell/face connectivity.
    - dgField stores one DG variable over the whole mesh, while each cell owns
      a GaussField that bundles interior and face Gauss-point values.
    - cellGaussField, faceGaussField, and boundaryGaussField behave much like
      scalar/vector/tensor containers, so many expressions can be written in a
      compact pointwise algebra style.
    - dgGeneralPDETerm assembles one semi-discrete DG operator for one cell.
    - dgTimeDiscretization drives the explicit stage loop and updates the
      registered conservative fields using the selected time scheme.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "dgGeomMesh.H"

// DG libs (must have)
#include "dgField.H"
#include "dgExpr.H"
#include "dgMath.H"
#include "dgCompressibleBoundaryManager.H"
#include "dgFluxSolverManager.H"
#include "dgLimiterManager.H"
#include "dgTimeDiscretization.H"
#include "dgThermoConservative.H"
#include "dgGeneralPDETerm.H"


// Test libs
#include "dgBasisField.H"
#include "debugProcessorFaceExchange.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    // Build the standard OpenFOAM run-time database:
    // - runTime reads system/controlDict,
    // - mesh reads the polyMesh and creates the underlying fvMesh.
    #include "createTime.H"
    #include "createMesh.H"

    // Build the DG-specific geometric layer on top of fvMesh. This is the
    // core geometry/data structure used throughout the solver: it provides
    // per-cell basis data, Gauss points, Jacobians, face connectivity, and
    // other DG geometric information needed for interpolation, integration,
    // flux assembly, and also the time-step estimate.
    #include "createDGMesh.H"

    // Create the DG solution fields and thermo model. The conservative
    // variables are the primary unknowns advanced by the time integrator,
    // while primitive/transport fields are reconstructed from them and used
    // to build fluxes and boundary states in a more natural form.
    #include "createDGFields.H"

    // Read dgSchemes and construct the two runtime "engines":
    // - dgFluxSolverManager selects numerical fluxes for each PDE term,
    // - dgTimeDiscretization manages stage looping and residual updates.
    #include "readSolverSettings.H"

    // Each outer iteration is one physical time step.
    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl
            << "deltaT = " << runTime.deltaTValue() << nl << endl;

        // The inner loop walks through explicit RK/Euler stages.
        // dgTimeDiscretization::loop() updates its internal stage index and
        // tells the solver when the current physical step is complete.
        while (timeDiscretization.loop())
        {
            // 1) Synchronize stage input fields and impose BCs.
            #include "prepareStageFields.H"
            // 2) Assemble the semi-discrete DG residual for this stage.
            //    This is where dgFoam feels most mathematical: Gauss fields
            //    support pointwise +, -, *, / just like ordinary values.
            #include "assembleStageResiduals.H"
            // 3) Apply the explicit stage update to the registered fields.
            #include "advanceStage.H"

            timeDiscretization.writeResiduals(Info);
        }

        // Clear stage bookkeeping so the next physical step starts from
        // stage 0 with fresh residual storage.
        timeDiscretization.reset();

        // Estimate the stable explicit time step from the current state
        // before incrementing the OpenFOAM time database.
        #include "calculateDeltaT.H"

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
