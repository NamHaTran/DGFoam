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
#include "PstreamReduceOps.H"
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


List<scalar> scalarResidualSummary(const List<List<scalar>>& residual)
{
    label maxDof = maxDofCount(residual);
    reduce(maxDof, maxOp<label>());

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
        reduce(maxResidual[dofI], maxOp<scalar>());
    }

    return maxResidual;
}


List<vector> vectorResidualSummary(const List<List<vector>>& residual)
{
    label maxDof = maxDofCount(residual);
    reduce(maxDof, maxOp<label>());

    List<vector> maxResidual(maxDof, Zero);

    forAll(residual, cellI)
    {
        const List<vector>& cellResidual = residual[cellI];

        forAll(cellResidual, dofI)
        {
            vector& maxResidualDof = maxResidual[dofI];
            const vector& residualDof = cellResidual[dofI];

            maxResidualDof.x() = max(maxResidualDof.x(), mag(residualDof.x()));
            maxResidualDof.y() = max(maxResidualDof.y(), mag(residualDof.y()));
            maxResidualDof.z() = max(maxResidualDof.z(), mag(residualDof.z()));
        }
    }

    forAll(maxResidual, dofI)
    {
        reduce(maxResidual[dofI].x(), maxOp<scalar>());
        reduce(maxResidual[dofI].y(), maxOp<scalar>());
        reduce(maxResidual[dofI].z(), maxOp<scalar>());
    }

    return maxResidual;
}


void accumulateScalarResidualReference
(
    List<scalar>& reference,
    const List<scalar>& summary
)
{
    if (!reference.size())
    {
        reference = summary;
        return;
    }

    if (reference.size() < summary.size())
    {
        const label oldSize = reference.size();
        reference.setSize(summary.size());

        for (label dofI = oldSize; dofI < summary.size(); ++dofI)
        {
            reference[dofI] = summary[dofI];
        }
    }

    forAll(summary, dofI)
    {
        reference[dofI] = max(reference[dofI], summary[dofI]);
    }
}


void accumulateVectorResidualReference
(
    List<scalar>& reference,
    const List<vector>& summary
)
{
    if (!reference.size())
    {
        reference.setSize(summary.size(), Zero);

        forAll(summary, dofI)
        {
            const vector& summaryDof = summary[dofI];
            reference[dofI] = max(summaryDof.x(), max(summaryDof.y(), summaryDof.z()));
        }

        return;
    }

    if (reference.size() < summary.size())
    {
        const label oldSize = reference.size();
        reference.setSize(summary.size(), Zero);

        for (label dofI = oldSize; dofI < summary.size(); ++dofI)
        {
            const vector& summaryDof = summary[dofI];
            reference[dofI] = max(summaryDof.x(), max(summaryDof.y(), summaryDof.z()));
        }
    }

    forAll(summary, dofI)
    {
        const vector& summaryDof = summary[dofI];
        const scalar summaryScale =
            max(summaryDof.x(), max(summaryDof.y(), summaryDof.z()));

        reference[dofI] = max(reference[dofI], summaryScale);
    }
}


scalar normalizeResidualValue(const scalar current, const scalar reference)
{
    return current/max(reference, VSMALL);
}


vector normalizeResidualValue(const vector& current, const scalar reference)
{
    return vector
    (
        normalizeResidualValue(current.x(), reference),
        normalizeResidualValue(current.y(), reference),
        normalizeResidualValue(current.z(), reference)
    );
}


void writeScalarResidualSummary
(
    Ostream& os,
    const word& fieldName,
    const List<List<scalar>>& residual,
    List<scalar>& reference,
    const bool referenceReady
)
{
    const List<scalar> maxResidual = scalarResidualSummary(residual);

    os  << "  scalar field " << fieldName << ':' << nl;

    if (!maxResidual.size())
    {
        os << "    no residual entries" << nl;
        return;
    }

    if (!referenceReady)
    {
        accumulateScalarResidualReference(reference, maxResidual);
    }

    const label dofI = 0;

    os << "    dof " << dofI << " : ";

    if (referenceReady)
    {
        os  << "max(|R|)/ref = "
            << normalizeResidualValue(maxResidual[dofI], reference[dofI]);
    }
    else
    {
        os << "max(|R|)/ref = 1";
    }

    os << nl;
}


void writeVectorResidualSummary
(
    Ostream& os,
    const word& fieldName,
    const List<List<vector>>& residual,
    List<scalar>& reference,
    const bool referenceReady
)
{
    const List<vector> maxResidual = vectorResidualSummary(residual);

    os  << "  vector field " << fieldName << ':' << nl;

    if (!maxResidual.size())
    {
        os << "    no residual entries" << nl;
        return;
    }

    if (!referenceReady)
    {
        accumulateVectorResidualReference(reference, maxResidual);
    }

    const label dofI = 0;

    os << "    dof " << dofI << " : ";

    if (referenceReady)
    {
        const vector scaledResidual =
            normalizeResidualValue(maxResidual[dofI], reference[dofI]);

        os  << "max(|Rx|)/ref = " << scaledResidual.x()
            << ", max(|Ry|)/ref = " << scaledResidual.y()
            << ", max(|Rz|)/ref = " << scaledResidual.z();
    }
    else
    {
        os  << "max(|Rx|)/ref = 1"
            << ", max(|Ry|)/ref = 1"
            << ", max(|Rz|)/ref = 1";
    }

    os << nl;
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
    timeStepStarted_(false),
    scalarResidualReference_(),
    vectorResidualReference_(),
    residualReferenceWarmupSteps_(5),
    completedTimeSteps_(0)
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
        << "  stages : " << nStage_ << nl << nl;
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

    const bool referenceReady =
        completedTimeSteps_ >= residualReferenceWarmupSteps_;

    const List<word>& scalarNames = timeScheme_().scalarFieldNames();
    const List<word>& vectorNames = timeScheme_().vectorFieldNames();

    if (scalarResidualReference_.size() != scalarNames.size())
    {
        scalarResidualReference_.setSize(scalarNames.size());
    }

    forAll(scalarResidualReference_, fieldI)
    {
        if (scalarResidualReference_[fieldI].size() != nStage_)
        {
            scalarResidualReference_[fieldI].setSize(nStage_);
        }
    }

    if (vectorResidualReference_.size() != vectorNames.size())
    {
        vectorResidualReference_.setSize(vectorNames.size());
    }

    forAll(vectorResidualReference_, fieldI)
    {
        if (vectorResidualReference_[fieldI].size() != nStage_)
        {
            vectorResidualReference_[fieldI].setSize(nStage_);
        }
    }

    forAll(scalarNames, fieldI)
    {
        writeScalarResidualSummary
        (
            os,
            scalarNames[fieldI],
            timeScheme_().scalarResidual(scalarNames[fieldI], stageI_),
            scalarResidualReference_[fieldI][stageI_],
            referenceReady
        );
    }

    forAll(vectorNames, fieldI)
    {
        writeVectorResidualSummary
        (
            os,
            vectorNames[fieldI],
            timeScheme_().vectorResidual(vectorNames[fieldI], stageI_),
            vectorResidualReference_[fieldI][stageI_],
            referenceReady
        );
    }
}

} // End namespace Foam

// ************************************************************************* //
