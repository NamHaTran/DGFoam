/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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

\*---------------------------------------------------------------------------*/

#include "dgTimeDiscretization.H"
#include "IOstreams.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

dgTimeDiscretization::dgTimeDiscretization
(
    const dictionary& dgSchemesDict,
    Time& runTime,
    const dgGeomMesh& mesh
)
:
    runTime_(runTime),
    mesh_(mesh),
    ddtSchemeDict_(dgSchemesDict.subDict("ddtSchemes")),
    schemeType_(word::null),
    timeScheme_(nullptr),
    nStage_(0),
    stageI_(0),
    timeStepStarted_(false)
{
    if (!ddtSchemeDict_.found("scheme"))
    {
        FatalIOErrorInFunction(dgSchemesDict)
            << "Missing entry 'ddtSchemes/scheme'." << nl
            << exit(FatalIOError);
    }

    schemeType_ = ddtSchemeDict_.get<word>("scheme");

    dictionary schemeDict(ddtSchemeDict_);

    // Optional coefficient block: ddtSchemeCoeffs/<SchemeName>Coeffs
    if (dgSchemesDict.found("ddtSchemeCoeffs"))
    {
        const dictionary& coeffsRoot = dgSchemesDict.subDict("ddtSchemeCoeffs");
        const word coeffKey(schemeType_ + "Coeffs");

        if (coeffsRoot.found(coeffKey))
        {
            schemeDict = coeffsRoot.subDict(coeffKey);
        }
    }

    timeScheme_.reset
    (
        dgTimeScheme::New("ddtScheme", schemeType_, schemeDict, mesh_).ptr()
    );

    // Cache stage count so the solver can drive the stage loop directly.
    nStage_ = timeScheme_().nStages();
    timeScheme_().setStage(stageI_);
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void dgTimeDiscretization::writeInfo(Ostream& os) const
{
    os  << "dgTimeDiscretization:" << nl
        << "  deltaT : " << runTime_.deltaTValue() << nl
        << "  scheme : " << schemeType_ << nl
        << "  stages : " << nStage_ << nl;
}

} // End namespace Foam

// ************************************************************************* //
