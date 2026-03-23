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

#include "troubleCellDetector.H"
#include "addToRunTimeSelectionTable.H"
#include "IOobject.H"
#include "dimensionedType.H"
#include "DynamicList.H"
#include "error.H"
#include "volFields.H"

namespace Foam
{

defineTypeNameAndDebug(troubleCellDetector, 0);
defineRunTimeSelectionTable(troubleCellDetector, dictionary);


word troubleCellDetector::lookupFieldName
(
    const dictionary& dict,
    const word& entryName,
    const word& defaultName
)
{
    if (dict.found(entryName))
    {
        return dict.get<word>(entryName);
    }

    return defaultName;
}


void troubleCellDetector::initializeCheckFields()
{
    const objectRegistry& obr = mesh_.getFvMesh();

    if (!dict_.found("checkFields"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required entry 'checkFields' in troubleCellDetector "
            << "dictionary." << nl
            << "Please provide a non-empty list of field names to be checked."
            << exit(FatalIOError);
    }

    const wordList entries = dict_.get<wordList>("checkFields");

    if (entries.empty())
    {
        FatalIOErrorInFunction(dict_)
            << "Entry 'checkFields' is empty in troubleCellDetector "
            << "dictionary." << nl
            << "Please provide at least one field name to be checked."
            << exit(FatalIOError);
    }

    DynamicList<checkedField> fields(entries.size());

    forAll(entries, i)
    {
        const word& entryName = entries[i];
        const word fieldName = lookupFieldName(dict_, entryName, entryName);

        checkedField fieldInfo;
        fieldInfo.entryName = entryName;
        fieldInfo.fieldName = fieldName;
        fieldInfo.isVector = false;
        fieldInfo.scalarFieldPtr = nullptr;
        fieldInfo.vectorFieldPtr = nullptr;

        if (obr.foundObject<dgField<scalar>>(fieldName))
        {
            fieldInfo.scalarFieldPtr =
                &obr.lookupObject<dgField<scalar>>(fieldName);
        }
        else if (obr.foundObject<dgField<vector>>(fieldName))
        {
            fieldInfo.isVector = true;
            fieldInfo.vectorFieldPtr =
                &obr.lookupObject<dgField<vector>>(fieldName);
        }
        else
        {
            FatalIOErrorInFunction(dict_)
                << "Field entry '" << entryName << "' resolved to object '"
                << fieldName << "', but no dgField<scalar> or dgField<vector>"
                << " with that name exists in the mesh registry."
                << exit(FatalIOError);
        }

        fields.append(fieldInfo);
    }

    checkFields_.transfer(fields);
}


troubleCellDetector::troubleCellDetector
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    report_(dict.lookupOrDefault<bool>("report", false)),
    checkFields_(),
    limitingIndicatorPtr_(nullptr)
{
    initializeCheckFields();

    if (report_)
    {
        limitingIndicatorPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "limitingIndicator",
                    mesh_.getFvMesh().time().timeName(),
                    mesh_.getFvMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_.getFvMesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        resetLimitingIndicator();
    }
}


void troubleCellDetector::setLimitingIndicator
(
    const label cellID,
    const bool flagged
) const
{
    if (!report_ || !limitingIndicatorPtr_.valid())
    {
        return;
    }

    limitingIndicatorPtr_().primitiveFieldRef()[cellID] =
        flagged ? scalar(1) : scalar(0);
}


void troubleCellDetector::resetLimitingIndicator() const
{
    if (!report_ || !limitingIndicatorPtr_.valid())
    {
        return;
    }

    auto& indicator = limitingIndicatorPtr_();
    indicator.primitiveFieldRef() = Zero;

    forAll(indicator.boundaryFieldRef(), patchI)
    {
        indicator.boundaryFieldRef()[patchI] = Zero;
    }
}


autoPtr<troubleCellDetector> troubleCellDetector::New
(
    const word& detectorType,
    const dictionary& dict,
    const dgGeomMesh& mesh
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->find(detectorType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown troubleCellDetector type: " << detectorType << nl
            << "Valid troubleCellDetector types are: "
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(dict, mesh);
}

} // End namespace Foam

// ************************************************************************* //
