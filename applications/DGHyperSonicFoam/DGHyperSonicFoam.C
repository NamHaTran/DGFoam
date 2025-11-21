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

// Test libs
#include "dgGeneralBoundaryManager.H"
#include "GaussField.H"
#include "dgField.H"
#include "dgBasisField.H"
#include "dgThermoConservative.H"

// math libs
#include "dgMath.H"


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

    // It's possible to iterate over every cell in a standard C++ for loop
    for (label cellI = 0; cellI < mesh.C().size(); cellI++)
    {
        Foam::dgBasisField basisField(cellI, dgMesh);

        if (cellI == 0)
        {
            const GaussField<scalar>& rhoG = rhoField.gaussFields()[cellI]; 
            const GaussField<vector>& rhoUG = rhoUField.gaussFields()[cellI];

            GaussField<vector>& UG = UField.gaussFields()[cellI];
            
            // Calculate U from rhoU and rho
            UG = rhoUG / rhoG;

            GaussField<scalar>& TG = thermo->T(cellI);
            GaussField<scalar>& CvG = thermo->Cv(cellI);
            GaussField<scalar>& muG = thermo->mu(cellI);

            Info << "U field : " << UG << endl;
            Info << "T field : " << TG << endl;
            Info << "Cv field : " << CvG << endl;
            Info << "mu field : " << muG << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
