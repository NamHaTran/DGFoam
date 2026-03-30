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

#include "positivityPreservingLimiter.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "dimensionedType.H"
#include "DynamicList.H"
#include "error.H"
#include "dgMath.H"
#include "volFields.H"

namespace Foam
{

defineTypeNameAndDebug(positivityPreservingLimiter, 0);
addToRunTimeSelectionTable(dgLimiter, positivityPreservingLimiter, dictionary);

namespace
{

void appendUniqueFieldName
(
    DynamicList<word>& fieldNames,
    const word& fieldName
)
{
    forAll(fieldNames, fieldI)
    {
        if (fieldNames[fieldI] == fieldName)
        {
            return;
        }
    }

    fieldNames.append(fieldName);
}


dictionary normalizePositivityLimiterDict(const dictionary& inputDict)
{
    if (!inputDict.found("density"))
    {
        FatalIOErrorInFunction(inputDict)
            << "Missing required entry 'density' in "
            << "positivityPreserving limiter dictionary."
            << exit(FatalIOError);
    }

    if (!inputDict.found("momentum"))
    {
        FatalIOErrorInFunction(inputDict)
            << "Missing required entry 'momentum' in "
            << "positivityPreserving limiter dictionary."
            << exit(FatalIOError);
    }

    if (!inputDict.found("energy"))
    {
        FatalIOErrorInFunction(inputDict)
            << "Missing required entry 'energy' in "
            << "positivityPreserving limiter dictionary."
            << exit(FatalIOError);
    }

    const word densityFieldName = inputDict.get<word>("density");
    const word momentumFieldName = inputDict.get<word>("momentum");
    const word energyFieldName = inputDict.get<word>("energy");

    if
    (
        densityFieldName == momentumFieldName
     || densityFieldName == energyFieldName
     || momentumFieldName == energyFieldName
    )
    {
        FatalIOErrorInFunction(inputDict)
            << "Entries 'density', 'momentum', and 'energy' must refer to "
            << "three distinct DG fields, but got density='"
            << densityFieldName << "', momentum='" << momentumFieldName
            << "', energy='" << energyFieldName << "'."
            << exit(FatalIOError);
    }

    DynamicList<word> mergedFieldNames(3);
    appendUniqueFieldName(mergedFieldNames, densityFieldName);
    appendUniqueFieldName(mergedFieldNames, momentumFieldName);
    appendUniqueFieldName(mergedFieldNames, energyFieldName);

    if (inputDict.found("fields"))
    {
        const wordList extraFieldNames = inputDict.get<wordList>("fields");

        forAll(extraFieldNames, fieldI)
        {
            appendUniqueFieldName(mergedFieldNames, extraFieldNames[fieldI]);
        }
    }

    dictionary normalizedDict(inputDict);
    normalizedDict.set("fields", wordList(mergedFieldNames));

    return normalizedDict;
}

} // End anonymous namespace


void positivityPreservingLimiter::initializeConservativeFields()
{
    densityFieldPtr_ = nullptr;
    momentumFieldPtr_ = nullptr;
    energyFieldPtr_ = nullptr;

    if (!dict_.found("density"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required entry 'density' in "
            << type() << " dictionary."
            << exit(FatalIOError);
    }

    if (!dict_.found("momentum"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required entry 'momentum' in "
            << type() << " dictionary."
            << exit(FatalIOError);
    }

    if (!dict_.found("energy"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required entry 'energy' in "
            << type() << " dictionary."
            << exit(FatalIOError);
    }

    densityFieldName_ = dict_.get<word>("density");
    momentumFieldName_ = dict_.get<word>("momentum");
    energyFieldName_ = dict_.get<word>("energy");

    if
    (
        densityFieldName_ == momentumFieldName_
     || densityFieldName_ == energyFieldName_
     || momentumFieldName_ == energyFieldName_
    )
    {
        FatalIOErrorInFunction(dict_)
            << "Entries 'density', 'momentum', and 'energy' must refer to "
            << "three distinct DG fields, but got density='"
            << densityFieldName_ << "', momentum='" << momentumFieldName_
            << "', energy='" << energyFieldName_ << "'."
            << exit(FatalIOError);
    }

    forAll(limitedFields_, fieldI)
    {
        const limitedField& fieldInfo = limitedFields_[fieldI];

        if (fieldInfo.fieldName == densityFieldName_)
        {
            if (fieldInfo.isVector)
            {
                FatalIOErrorInFunction(dict_)
                    << "Density field '" << densityFieldName_
                    << "' must be a dgField<scalar>."
                    << exit(FatalIOError);
            }

            densityFieldPtr_ = fieldInfo.scalarFieldPtr;
            continue;
        }

        if (fieldInfo.fieldName == momentumFieldName_)
        {
            if (!fieldInfo.isVector)
            {
                FatalIOErrorInFunction(dict_)
                    << "Momentum field '" << momentumFieldName_
                    << "' must be a dgField<vector>."
                    << exit(FatalIOError);
            }

            momentumFieldPtr_ = fieldInfo.vectorFieldPtr;
            continue;
        }

        if (fieldInfo.fieldName == energyFieldName_)
        {
            if (fieldInfo.isVector)
            {
                FatalIOErrorInFunction(dict_)
                    << "Energy field '" << energyFieldName_
                    << "' must be a dgField<scalar>."
                    << exit(FatalIOError);
            }

            energyFieldPtr_ = fieldInfo.scalarFieldPtr;
        }
    }

    if (!densityFieldPtr_)
    {
        FatalIOErrorInFunction(dict_)
            << "Density field '" << densityFieldName_
            << "' could not be resolved from the limiter field list."
            << exit(FatalIOError);
    }

    if (!momentumFieldPtr_)
    {
        FatalIOErrorInFunction(dict_)
            << "Momentum field '" << momentumFieldName_
            << "' could not be resolved from the limiter field list."
            << exit(FatalIOError);
    }

    if (!energyFieldPtr_)
    {
        FatalIOErrorInFunction(dict_)
            << "Energy field '" << energyFieldName_
            << "' could not be resolved from the limiter field list."
            << exit(FatalIOError);
    }
}


scalar positivityPreservingLimiter::densityMean(const label cellID) const
{
    if (!densityFieldPtr_->hasDof())
    {
        FatalErrorInFunction
            << type() << " requires a dofField for density field '"
            << densityFieldPtr_->name() << "'."
            << abort(FatalError);
    }

    const List<scalar>& rhoDof = densityFieldPtr_->dof()[cellID].dof();

    if (rhoDof.empty())
    {
        FatalErrorInFunction
            << "Density field '" << densityFieldPtr_->name()
            << "' has no modal coefficients in cell " << cellID << '.'
            << abort(FatalError);
    }

    return rhoDof[0];
}


scalar positivityPreservingLimiter::modifiedDensityValue
(
    const scalar rhoValue,
    const scalar meanRho,
    const scalar theta1
) const
{
    return meanRho + theta1*(rhoValue - meanRho);
}


scalar positivityPreservingLimiter::calcRhoThermoEnergy
(
    const scalar rho,
    const vector& rhoU,
    const scalar E
) const
{
    const scalar rhoSafe =
        (mag(rho) > epsilon_ ? rho : (rho >= 0 ? epsilon_ : -epsilon_));

    return E - 0.5*magSqr(rhoU)/rhoSafe;
}


scalar positivityPreservingLimiter::thetaCoeff
(
    const scalar meanValue,
    const scalar minValue,
    const scalar omega
) const
{
    const scalar scale = max(scalar(1), mag(meanValue));

    if (mag(meanValue - minValue) <= epsilon_*scale)
    {
        return 1.0;
    }

    const scalar rawTheta = (meanValue - omega)/(meanValue - minValue);
    return min(max(rawTheta, scalar(0)), scalar(1));
}


void positivityPreservingLimiter::computeThetaCoeffs
(
    const label cellID,
    scalar& theta1,
    scalar& theta2
) const
{
    const scalar meanRho = densityMean(cellID);
    const dgGeomCell& cell = *mesh_.cells()[cellID];

    const GaussField<scalar>& rhoGF = densityFieldPtr_->gaussFields()[cellID];
    const GaussField<vector>& rhoUGF = momentumFieldPtr_->gaussFields()[cellID];
    const GaussField<scalar>& EGF = energyFieldPtr_->gaussFields()[cellID];

    const cellGaussField<scalar>& rhoCell = rhoGF.cellField();
    const cellGaussField<vector>& rhoUCell = rhoUGF.cellField();
    const cellGaussField<scalar>& ECell = EGF.cellField();

    const faceGaussField<scalar>& rhoFace = rhoGF.faceField();
    const faceGaussField<vector>& rhoUFace = rhoUGF.faceField();
    const faceGaussField<scalar>& EFace = EGF.faceField();

    const List<scalar>& cellWeights = cell.weights();
    const List<scalar>& J3D = cell.J3D();

    label nFaceSamples = 0;

    for (label localFaceI = 0; localFaceI < rhoFace.nFaces(); ++localFaceI)
    {
        nFaceSamples += rhoFace.faces()[localFaceI]->weights().size();
    }

    DynamicList<conservativeSample> samples(rhoCell.size() + nFaceSamples);
    scalar minRho = GREAT;

    forAll(rhoCell, gpI)
    {
        conservativeSample sample;
        sample.rho = rhoCell[gpI];
        sample.rhoU = rhoUCell[gpI];
        sample.E = ECell[gpI];
        sample.meanWeight = cellWeights[gpI]*J3D[gpI];
        sample.contributesToMean = true;

        minRho = min(minRho, sample.rho);
        samples.append(sample);
    }

    for (label localFaceI = 0; localFaceI < rhoFace.nFaces(); ++localFaceI)
    {
        const dgGeomFace& face = *rhoFace.faces()[localFaceI];
        const List<scalar>& weights = face.weights();

        forAll(weights, gpI)
        {
            conservativeSample sample;
            sample.rho = rhoFace.minusValueOnFace(localFaceI, gpI);
            sample.rhoU = rhoUFace.minusValueOnFace(localFaceI, gpI);
            sample.E = EFace.minusValueOnFace(localFaceI, gpI);
            sample.meanWeight = Zero;
            sample.contributesToMean = false;

            minRho = min(minRho, sample.rho);
            samples.append(sample);
        }
    }

    theta1 = thetaCoeff(meanRho, minRho, min(epsilon_, meanRho));

    scalar weightedRhoThermoEnergySum = Zero;
    scalar totalVolume = Zero;
    scalar minRhoThermoEnergy = GREAT;

    forAll(samples, sampleI)
    {
        const conservativeSample& sample = samples[sampleI];
        const scalar rho =
            modifiedDensityValue(sample.rho, meanRho, theta1);
        const scalar rhoThermoEnergy =
            calcRhoThermoEnergy(rho, sample.rhoU, sample.E);

        minRhoThermoEnergy =
            min(minRhoThermoEnergy, rhoThermoEnergy);

        if (sample.contributesToMean)
        {
            weightedRhoThermoEnergySum +=
                sample.meanWeight*rhoThermoEnergy;
            totalVolume += sample.meanWeight;
        }
    }

    const scalar meanRhoThermoEnergy =
        weightedRhoThermoEnergySum/max(totalVolume, VSMALL);

    theta2 =
        thetaCoeff
        (
            meanRhoThermoEnergy,
            minRhoThermoEnergy,
            min(min(epsilon_, meanRho), meanRhoThermoEnergy)
        );
}


void positivityPreservingLimiter::scaleHighOrderModes
(
    dgField<scalar>& field,
    const label cellID,
    const scalar theta
) const
{
    if (!field.hasDof() || theta >= 1.0)
    {
        return;
    }

    cellDof<scalar>& cellModes = field.dof()[cellID];

    for (label modeI = 1; modeI < cellModes.nDof(); ++modeI)
    {
        cellModes[modeI] *= theta;
    }

    field.dof().updateCellDof(cellID);
}


void positivityPreservingLimiter::scaleHighOrderModes
(
    dgField<vector>& field,
    const label cellID,
    const scalar theta
) const
{
    if (!field.hasDof() || theta >= 1.0)
    {
        return;
    }

    cellDof<vector>& cellModes = field.dof()[cellID];

    for (label modeI = 1; modeI < cellModes.nDof(); ++modeI)
    {
        cellModes[modeI] *= theta;
    }

    field.dof().updateCellDof(cellID);
}


void positivityPreservingLimiter::resetThetaFields() const
{
    if (!writeTheta_)
    {
        return;
    }

    if (theta1FieldPtr_.valid())
    {
        auto& theta1Field = theta1FieldPtr_();
        theta1Field.primitiveFieldRef() = scalar(1);

        forAll(theta1Field.boundaryFieldRef(), patchI)
        {
            theta1Field.boundaryFieldRef()[patchI] = scalar(1);
        }
    }

    if (theta2FieldPtr_.valid())
    {
        auto& theta2Field = theta2FieldPtr_();
        theta2Field.primitiveFieldRef() = scalar(1);

        forAll(theta2Field.boundaryFieldRef(), patchI)
        {
            theta2Field.boundaryFieldRef()[patchI] = scalar(1);
        }
    }
}


void positivityPreservingLimiter::setThetaFields
(
    const label cellID,
    const scalar theta1,
    const scalar theta2
) const
{
    if (!writeTheta_)
    {
        return;
    }

    if (theta1FieldPtr_.valid())
    {
        theta1FieldPtr_().primitiveFieldRef()[cellID] = theta1;
    }

    if (theta2FieldPtr_.valid())
    {
        theta2FieldPtr_().primitiveFieldRef()[cellID] = theta2;
    }
}


void positivityPreservingLimiter::postCorrect()
{
    dgLimiter::postCorrect();

    label nLimitedCells = nLimitedCells_;
    reduce(nLimitedCells, sumOp<label>());

    Info<< "Limiter " << type()
        << " applied on " << nLimitedCells
        << " cell(s)." << nl;
}


void positivityPreservingLimiter::limitField
(
    dgField<scalar>& field,
    const label cellID
)
{
    const scalar theta =
        (&field == densityFieldPtr_ ? theta1_*theta2_ : theta2_);

    scaleHighOrderModes(field, cellID, theta);
}


void positivityPreservingLimiter::limitField
(
    dgField<vector>& field,
    const label cellID
)
{
    scaleHighOrderModes(field, cellID, theta2_);
}


positivityPreservingLimiter::positivityPreservingLimiter
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgLimiter(normalizePositivityLimiterDict(dict), mesh),
    densityFieldName_(word::null),
    momentumFieldName_(word::null),
    energyFieldName_(word::null),
    densityFieldPtr_(nullptr),
    momentumFieldPtr_(nullptr),
    energyFieldPtr_(nullptr),
    epsilon_(dict.lookupOrDefault<scalar>("tolerance", 1.0e-8)),
    theta1_(1.0),
    theta2_(1.0),
    writeTheta_(dict.lookupOrDefault<bool>("writeTheta", false)),
    theta1FieldPtr_(nullptr),
    theta2FieldPtr_(nullptr)
{
    initializeConservativeFields();

    if (writeTheta_)
    {
        theta1FieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "theta1",
                    mesh_.getFvMesh().time().timeName(),
                    mesh_.getFvMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_.getFvMesh(),
                dimensionedScalar("one", dimless, scalar(1))
            )
        );

        theta2FieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "theta2",
                    mesh_.getFvMesh().time().timeName(),
                    mesh_.getFvMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_.getFvMesh(),
                dimensionedScalar("one", dimless, scalar(1))
            )
        );

        resetThetaFields();
    }
}


void positivityPreservingLimiter::read(const dictionary& dict)
{
    dgLimiter::read(normalizePositivityLimiterDict(dict));
    epsilon_ = dict.lookupOrDefault<scalar>("tolerance", 1.0e-8);
    theta1_ = 1.0;
    theta2_ = 1.0;
    writeTheta_ = dict.lookupOrDefault<bool>("writeTheta", false);
    theta1FieldPtr_.clear();
    theta2FieldPtr_.clear();
    initializeConservativeFields();

    if (writeTheta_)
    {
        theta1FieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "theta1",
                    mesh_.getFvMesh().time().timeName(),
                    mesh_.getFvMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_.getFvMesh(),
                dimensionedScalar("one", dimless, scalar(1))
            )
        );

        theta2FieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "theta2",
                    mesh_.getFvMesh().time().timeName(),
                    mesh_.getFvMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_.getFvMesh(),
                dimensionedScalar("one", dimless, scalar(1))
            )
        );

        resetThetaFields();
    }
}


void positivityPreservingLimiter::correct()
{
    preCorrect();
    resetThetaFields();

    nLimitedCells_ = 0;

    if (hasDetector())
    {
        detector().resetLimitingIndicator();
    }

    for (label cellID = 0; cellID < mesh_.nCells(); ++cellID)
    {
        computeThetaCoeffs(cellID, theta1_, theta2_);
        setThetaFields(cellID, theta1_, theta2_);

        const scalar densityTheta = theta1_*theta2_;
        const scalar otherTheta = theta2_;

        if (densityTheta >= 1.0 && otherTheta >= 1.0)
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
