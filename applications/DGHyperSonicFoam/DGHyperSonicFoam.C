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
#include "dgRefFace.H"
#include "dgRefCell.H"
#include "dgGeomCell.H"
#include "dgGeomMesh.H"

// Test libs
#include "dgGeneralBoundaryManager.H"
#include "GaussField.H"
#include "dgField.H"
#include "dgBasisField.H"

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

        // Thermodynamic properties:
        const scalar R = 287.0;          // Specific gas constant for air
        const scalar gamma = 1.4;        // Specific heat ratio for air


        if (cellI == 0)
        {
            // Test lookup dgField and use
            const dgField<vector>& UDG = dgMesh.getFvMesh().lookupObject<dgField<vector>>("U");

            GaussField<vector> UGauss = UDG.gaussFields()[cellI];

            Info << "UGauss at cell " << cellI << ": " << UGauss << endl;

            GaussField<scalar> pGauss = pField.gaussFields()[cellI];
            GaussField<scalar> TGauss = TField.gaussFields()[cellI];

            GaussField<scalar> rhoGauss = pGauss / (R * TGauss);
            GaussField<scalar> eGauss = TGauss * gamma / (gamma - 1.0);
            GaussField<scalar> aGauss = sqrt(gamma * pGauss / rhoGauss);
            GaussField<scalar> machGauss = mag(UGauss) / aGauss;
            GaussField<scalar> EGauss = eGauss + 0.5 * magSqr(UGauss);

            Info << "rhoGauss at cell " << cellI << ": " << rhoGauss << endl;
            Info << "eGauss at cell " << cellI << ": " << eGauss << endl;
            Info << "aGauss at cell " << cellI << ": " << aGauss << endl;
            Info << "machGauss at cell " << cellI << ": " << machGauss << endl;
            Info << "EGauss at cell " << cellI << ": " << EGauss << endl;

            /*
            Foam::GaussField<tensor> T1
            (
                cellI,
                &dgMesh,
                Foam::tensor
                ( 1, 2, 3, 
                4, 4, 5,
                7, 3, 1 )
            );
            
            Foam::GaussField<vector> V1
            (
                cellI,
                &dgMesh,
                Foam::vector(1, 2, 3)
            );

            Foam::GaussField<symmTensor> symmT1
            (
                cellI,
                &dgMesh,
                Foam::symmTensor(1.0, 0.2, 0.3, 2.0, 0.4, 3.0)
            );

            Foam::GaussField<symmTensor> symmT2 = symmT1*2.0;

            Foam::GaussField<tensor> T2
            (
                cellI,
                &dgMesh,
                Foam::tensor

                ( 1.0, 0.2, 0.3,
                0.2, 2.0, 0.4,
                0.3, 0.4, 3.0 )
            );
            */
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
