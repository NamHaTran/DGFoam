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

namespace
{

label maxDofCount(const List<List<scalar>>& residual)
{
    label maxDof = 0;

    forAll(residual, cellI)
    {
        maxDof = max(maxDof, residual[cellI].size());
    }

    return maxDof;
}


label maxDofCount(const List<List<vector>>& residual)
{
    label maxDof = 0;

    forAll(residual, cellI)
    {
        maxDof = max(maxDof, residual[cellI].size());
    }

    return maxDof;
}


void writeScalarResidualSummary
(
    Ostream& os,
    const word& fieldName,
    const List<List<scalar>>& residual
)
{
    const label maxDof = maxDofCount(residual);

    os  << "  scalar field " << fieldName << ':' << nl;

    if (!maxDof)
    {
        os << "    no residual entries" << nl;
        return;
    }

    List<scalar> maxResidual(maxDof, Zero);

    forAll(residual, cellI)
    {
        const List<scalar>& cellResidual = residual[cellI];

        forAll(cellResidual, dofI)
        {
            maxResidual[dofI] = max
            (
                maxResidual[dofI],
                mag(cellResidual[dofI])
            );
        }
    }

    forAll(maxResidual, dofI)
    {
        os  << "    dof " << dofI
            << " : max(|R|) = " << maxResidual[dofI] << nl;
    }
}


void writeVectorResidualSummary
(
    Ostream& os,
    const word& fieldName,
    const List<List<vector>>& residual
)
{
    const label maxDof = maxDofCount(residual);

    os  << "  vector field " << fieldName << ':' << nl;

    if (!maxDof)
    {
        os << "    no residual entries" << nl;
        return;
    }

    List<scalar> maxResidualX(maxDof, Zero);
    List<scalar> maxResidualY(maxDof, Zero);
    List<scalar> maxResidualZ(maxDof, Zero);

    forAll(residual, cellI)
    {
        const List<vector>& cellResidual = residual[cellI];

        forAll(cellResidual, dofI)
        {
            const vector& residualDof = cellResidual[dofI];

            maxResidualX[dofI] = max(maxResidualX[dofI], mag(residualDof.x()));
            maxResidualY[dofI] = max(maxResidualY[dofI], mag(residualDof.y()));
            maxResidualZ[dofI] = max(maxResidualZ[dofI], mag(residualDof.z()));
        }
    }

    forAll(maxResidualX, dofI)
    {
        os  << "    dof " << dofI
            << " : max(|Rx|) = " << maxResidualX[dofI]
            << ", max(|Ry|) = " << maxResidualY[dofI]
            << ", max(|Rz|) = " << maxResidualZ[dofI] << nl;
    }
}

} // End anonymous namespace

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
        << "  scheme : " << schemeType_ << nl
        << "  stages : " << nStage_ << nl;
}


void dgTimeDiscretization::writeResiduals(Ostream& os) const
{
    if (stageI_ < 0 || stageI_ >= nStage_)
    {
        FatalErrorInFunction
            << "Stage index " << stageI_
            << " is out of range [0," << (nStage_ - 1) << "]."
            << exit(FatalError);
    }

    os  << "dgTimeDiscretization residuals:" << nl
        << "  stage : " << stageI_ << nl;

    const List<word>& scalarNames = timeScheme_().scalarFieldNames();
    const List<word>& vectorNames = timeScheme_().vectorFieldNames();

    forAll(scalarNames, fieldI)
    {
        writeScalarResidualSummary
        (
            os,
            scalarNames[fieldI],
            timeScheme_().scalarResidual(scalarNames[fieldI], stageI_)
        );
    }

    forAll(vectorNames, fieldI)
    {
        writeVectorResidualSummary
        (
            os,
            vectorNames[fieldI],
            timeScheme_().vectorResidual(vectorNames[fieldI], stageI_)
        );
    }
}

} // End namespace Foam

// ************************************************************************* //
