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

#include "dgLimiter.H"
#include "addToRunTimeSelectionTable.H"
#include "IOobject.H"
#include "dimensionedType.H"
#include "DynamicList.H"
#include "error.H"
#include "volFields.H"

namespace Foam
{

defineTypeNameAndDebug(dgLimiter, 0);
defineRunTimeSelectionTable(dgLimiter, dictionary);


void dgLimiter::initializeLimitedFields()
{
    const objectRegistry& obr = mesh_.getFvMesh();

    if (!dict_.found("fields"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required entry 'fields' in dgLimiter dictionary." << nl
            << "Please provide a non-empty list of DG field names to limit."
            << exit(FatalIOError);
    }

    const wordList fieldNames = dict_.get<wordList>("fields");

    if (fieldNames.empty())
    {
        FatalIOErrorInFunction(dict_)
            << "Entry 'fields' is empty in dgLimiter dictionary." << nl
            << "Please provide at least one DG field name to limit."
            << exit(FatalIOError);
    }

    DynamicList<limitedField> fields(fieldNames.size());

    forAll(fieldNames, fieldI)
    {
        const word& fieldName = fieldNames[fieldI];

        limitedField fieldInfo;
        fieldInfo.fieldName = fieldName;
        fieldInfo.isVector = false;
        fieldInfo.scalarFieldPtr = nullptr;
        fieldInfo.vectorFieldPtr = nullptr;

        if (obr.foundObject<dgField<scalar>>(fieldName))
        {
            fieldInfo.scalarFieldPtr =
                &obr.lookupObjectRef<dgField<scalar>>(fieldName);
        }
        else if (obr.foundObject<dgField<vector>>(fieldName))
        {
            fieldInfo.isVector = true;
            fieldInfo.vectorFieldPtr =
                &obr.lookupObjectRef<dgField<vector>>(fieldName);
        }
        else
        {
            FatalIOErrorInFunction(dict_)
                << "Field '" << fieldName
                << "' does not exist in the mesh registry as either "
                << "dgField<scalar> or dgField<vector>."
                << exit(FatalIOError);
        }

        fields.append(fieldInfo);
    }

    limitedFields_.transfer(fields);
}


void dgLimiter::initializeDetector()
{
    if (!dict_.found("detector"))
    {
        detector_.clear();
        return;
    }

    const word detectorType = dict_.get<word>("detector");
    const word detectorCoeffsName = detectorType + "Coeffs";

    if (!dict_.found(detectorCoeffsName))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required sub-dictionary '" << detectorCoeffsName
            << "' for detector type '" << detectorType << "'."
            << exit(FatalIOError);
    }

    detector_ =
        troubleCellDetector::New
        (
            detectorType,
            dict_.subDict(detectorCoeffsName),
            mesh_
        );
}


void dgLimiter::limitCell
(
    const limitedField& fieldInfo,
    const label cellID
)
{
    if (fieldInfo.isVector)
    {
        limitField(*fieldInfo.vectorFieldPtr, cellID);
        return;
    }

    limitField(*fieldInfo.scalarFieldPtr, cellID);
}


void dgLimiter::preCorrect()
{
    limitedCellIDs_.clear();
    resetLimitedCellIndicator();
}


void dgLimiter::postCorrect()
{
    rebuildLimitedCellGaussFields();
}


void dgLimiter::cacheLimitedCell(const label cellID)
{
    limitedCellIDs_.append(cellID);
    setLimitedCellIndicator(cellID, true);
}


void dgLimiter::setLimitedCellIndicator
(
    const label cellID,
    const bool limited
) const
{
    if (!report_ || !limitedCellFieldPtr_.valid())
    {
        return;
    }

    limitedCellFieldPtr_().primitiveFieldRef()[cellID] =
        limited ? scalar(1) : scalar(0);
}


void dgLimiter::resetLimitedCellIndicator() const
{
    if (!report_ || !limitedCellFieldPtr_.valid())
    {
        return;
    }

    auto& indicator = limitedCellFieldPtr_();
    indicator.primitiveFieldRef() = Zero;

    forAll(indicator.boundaryFieldRef(), patchI)
    {
        indicator.boundaryFieldRef()[patchI] = Zero;
    }
}


void dgLimiter::rebuildLimitedCellGaussFields()
{
    if (!limitedCellIDs_.size())
    {
        return;
    }

    List<bool> rebuildMask(mesh_.nCells(), false);
    DynamicList<label> rebuildCellIDs(limitedCellIDs_.size()*2);

    forAll(limitedCellIDs_, limitedCellI)
    {
        const label cellID = limitedCellIDs_[limitedCellI];

        if (!rebuildMask[cellID])
        {
            rebuildMask[cellID] = true;
            rebuildCellIDs.append(cellID);
        }

        const labelList& neighbours = mesh_.cells()[cellID]->neighborCells();

        forAll(neighbours, neighI)
        {
            const label neighID = neighbours[neighI];

            if (neighID < 0 || rebuildMask[neighID])
            {
                continue;
            }

            rebuildMask[neighID] = true;
            rebuildCellIDs.append(neighID);
        }
    }

    forAll(limitedFields_, fieldI)
    {
        const limitedField& fieldInfo = limitedFields_[fieldI];

        if (fieldInfo.isVector)
        {
            dgField<vector>& field = *fieldInfo.vectorFieldPtr;

            if (!field.hasDof())
            {
                continue;
            }

            forAll(rebuildCellIDs, rebuildCellI)
            {
                const label cellID = rebuildCellIDs[rebuildCellI];
                field.gaussFields()[cellID].interpolateFromDof();
            }
        }
        else
        {
            dgField<scalar>& field = *fieldInfo.scalarFieldPtr;

            if (!field.hasDof())
            {
                continue;
            }

            forAll(rebuildCellIDs, rebuildCellI)
            {
                const label cellID = rebuildCellIDs[rebuildCellI];
                field.gaussFields()[cellID].interpolateFromDof();
            }
        }
    }
}


dgLimiter::dgLimiter
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    limitedFields_(),
    detector_(nullptr),
    dict_(dict),
    mesh_(mesh),
    nLimitedCells_(0),
    report_(dict.lookupOrDefault<bool>("report", false)),
    limitedCellFieldPtr_(nullptr)
{
    initializeLimitedFields();
    initializeDetector();

    if (report_)
    {
        limitedCellFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "limitedCell",
                    mesh_.getFvMesh().time().timeName(),
                    mesh_.getFvMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_.getFvMesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        resetLimitedCellIndicator();
    }
}


autoPtr<dgLimiter> dgLimiter::New
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
{
    if (!dict.found("type"))
    {
        FatalIOErrorInFunction(dict)
            << "Missing required entry 'type' in dgLimiter dictionary."
            << exit(FatalIOError);
    }

    const word limiterType = dict.get<word>("type");
    auto cstrIter = dictionaryConstructorTablePtr_->find(limiterType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown dgLimiter type: " << limiterType << nl
            << "Valid dgLimiter types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(dict, mesh);
}


void dgLimiter::read(const dictionary& dict)
{
    dict_ = dict;
    limitedFields_.clear();
    detector_.clear();
    nLimitedCells_ = 0;
    report_ = dict.lookupOrDefault<bool>("report", false);
    limitedCellFieldPtr_.clear();

    initializeLimitedFields();
    initializeDetector();

    if (report_)
    {
        limitedCellFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "limitedCell",
                    mesh_.getFvMesh().time().timeName(),
                    mesh_.getFvMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_.getFvMesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        resetLimitedCellIndicator();
    }
}


void dgLimiter::correct()
{
    if (!detector_.valid())
    {
        FatalErrorInFunction
            << "Limiter '" << type()
            << "' requires a configured trouble-cell detector before "
            << "calling correct()."
            << abort(FatalError);
    }

    preCorrect();

    nLimitedCells_ = 0;
    detector_->resetLimitingIndicator();

    for (label cellID = 0; cellID < mesh_.nCells(); ++cellID)
    {
        if (!detector_->detect(cellID))
        {
            continue;
        }

        ++nLimitedCells_;
        cacheLimitedCell(cellID);

        forAll(limitedFields_, fieldI)
        {
            limitCell(limitedFields_[fieldI], cellID);
        }
    }

    postCorrect();
}

} // End namespace Foam

// ************************************************************************* //
