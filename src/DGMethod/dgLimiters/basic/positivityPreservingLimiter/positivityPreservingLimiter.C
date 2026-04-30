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


void positivityPreservingLimiter::initializeThermo()
{
    thermoPtr_ = nullptr;

    const objectRegistry& obr = mesh_.getFvMesh();

    if (!obr.foundObject<dgThermoConservative>("dgThermoConservative"))
    {
        FatalErrorInFunction
            << "Could not find dgThermoConservative in the mesh registry while "
            << "constructing limiter '" << type() << "'."
            << abort(FatalError);
    }

    thermoPtr_ = &obr.lookupObject<dgThermoConservative>("dgThermoConservative");
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


scalar positivityPreservingLimiter::calcRhoThermoEnergyFloor
(
    const scalar rho
) const
{
    const scalar rhoFloor = max(rho, rhoMin_);
    const scalar heFloor = thermoPtr_->calcHeFromRhoT(rhoFloor, Tmin_);
    scalar floor = rhoFloor*heFloor;

    if (thermoPtr_->heIsEnthalpy())
    {
        floor -= thermoPtr_->eos().calcPFromRhoT(rhoFloor, Tmin_);
    }

    return floor;
}


scalar positivityPreservingLimiter::safeDensityForDivision
(
    const scalar rho
) const
{
    return
        mag(rho) > epsilon_
      ? rho
      : (rho >= 0 ? epsilon_ : -epsilon_);
}


scalar positivityPreservingLimiter::calcMeanTemperature
(
    const label cellID,
    const scalar rho,
    const vector& rhoU,
    const scalar E
) const
{
    const scalar rhoSafe = safeDensityForDivision(rho);
    const vector U = rhoU/rhoSafe;
    const scalar he = E/rhoSafe - 0.5*magSqr(U);

    return thermoPtr_->calcTemperatureFromRhoHe(cellID, rhoSafe, he);
}


scalar positivityPreservingLimiter::calcMeanTotalEnergy
(
    const scalar rho,
    const vector& U,
    const scalar T
) const
{
    const scalar he = thermoPtr_->calcHeFromRhoT(rho, T);
    scalar E = rho*(he + 0.5*magSqr(U));

    if (thermoPtr_->heIsEnthalpy())
    {
        E -= thermoPtr_->eos().calcPFromRhoT(rho, T);
    }

    return E;
}


bool positivityPreservingLimiter::limitCellMeanValue(const label cellID)
{
    if
    (
        !densityFieldPtr_->hasDof()
     || !momentumFieldPtr_->hasDof()
     || !energyFieldPtr_->hasDof()
    )
    {
        return false;
    }

    cellDof<scalar>& rhoModes = densityFieldPtr_->dof()[cellID];
    cellDof<vector>& rhoUModes = momentumFieldPtr_->dof()[cellID];
    cellDof<scalar>& EModes = energyFieldPtr_->dof()[cellID];

    if
    (
        rhoModes.nDof() < 1
     || rhoUModes.nDof() < 1
     || EModes.nDof() < 1
    )
    {
        return false;
    }

    bool limited = false;

    const scalar rho0NoClamp = rhoModes[0];
    const scalar rho0Clamped = min(max(rho0NoClamp, rhoMin_), rhoMax_);

    if (rho0Clamped != rho0NoClamp)
    {
        rhoModes[0] = rho0Clamped;
        limited = true;
    }

    const vector U0 = rhoUModes[0]/safeDensityForDivision(rho0NoClamp);
    const scalar magU0 = mag(U0);

    if (magU0 > magUMax_ && magU0 > SMALL)
    {
        rhoUModes[0] = rhoModes[0]*U0*(magUMax_/magU0);
        limited = true;
    }

    const scalar T0 =
        calcMeanTemperature(cellID, rhoModes[0], rhoUModes[0], EModes[0]);
    const scalar TClamped = min(max(T0, Tmin_), Tmax_);

    if (TClamped != T0)
    {
        const vector UMean =
            rhoUModes[0]/safeDensityForDivision(rhoModes[0]);
        EModes[0] = calcMeanTotalEnergy(rhoModes[0], UMean, TClamped);
        limited = true;
    }

    return limited;
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

    const GaussField<scalar>& rhoGF = densityFieldPtr_->gaussFields()[cellID];
    const GaussField<vector>& rhoUGF = momentumFieldPtr_->gaussFields()[cellID];
    const GaussField<scalar>& EGF = energyFieldPtr_->gaussFields()[cellID];

    const cellGaussField<scalar>& rhoCell = rhoGF.cellField();
    const cellGaussField<vector>& rhoUCell = rhoUGF.cellField();
    const cellGaussField<scalar>& ECell = EGF.cellField();

    const faceGaussField<scalar>& rhoFace = rhoGF.faceField();
    const faceGaussField<vector>& rhoUFace = rhoUGF.faceField();
    const faceGaussField<scalar>& EFace = EGF.faceField();

    scalar minRho = GREAT;

    forAll(rhoCell, gpI)
    {
        minRho = min(minRho, rhoCell[gpI]);
    }

    for (label localFaceI = 0; localFaceI < rhoFace.nFaces(); ++localFaceI)
    {
        const dgGeomFace& face = *rhoFace.faces()[localFaceI];
        const List<scalar>& weights = face.weights();

        forAll(weights, gpI)
        {
            minRho =
                min(minRho, rhoFace.minusValueOnFace(localFaceI, gpI));
        }
    }

    theta1 = thetaCoeff(meanRho, minRho, min(rhoMin_, meanRho));

    scalar minRhoThermoEnergyMargin = GREAT;

    forAll(rhoCell, gpI)
    {
        const scalar rho =
            modifiedDensityValue(rhoCell[gpI], meanRho, theta1);
        const scalar rhoThermoEnergy =
            calcRhoThermoEnergy(rho, rhoUCell[gpI], ECell[gpI]);
        const scalar rhoThermoEnergyFloor =
            calcRhoThermoEnergyFloor(rho);

        minRhoThermoEnergyMargin =
            min
            (
                minRhoThermoEnergyMargin,
                rhoThermoEnergy - rhoThermoEnergyFloor
            );
    }

    for (label localFaceI = 0; localFaceI < rhoFace.nFaces(); ++localFaceI)
    {
        const dgGeomFace& face = *rhoFace.faces()[localFaceI];
        const List<scalar>& weights = face.weights();

        forAll(weights, gpI)
        {
            const scalar rho =
                modifiedDensityValue
                (
                    rhoFace.minusValueOnFace(localFaceI, gpI),
                    meanRho,
                    theta1
                );
            const scalar rhoThermoEnergy =
                calcRhoThermoEnergy
                (
                    rho,
                    rhoUFace.minusValueOnFace(localFaceI, gpI),
                    EFace.minusValueOnFace(localFaceI, gpI)
                );
            const scalar rhoThermoEnergyFloor =
                calcRhoThermoEnergyFloor(rho);

            minRhoThermoEnergyMargin =
                min
                (
                    minRhoThermoEnergyMargin,
                    rhoThermoEnergy - rhoThermoEnergyFloor
                );
        }
    }

    const vector& meanRhoU = momentumFieldPtr_->dof()[cellID].dof()[0];
    const scalar meanE = energyFieldPtr_->dof()[cellID].dof()[0];
    const scalar meanRhoThermoEnergy =
        calcRhoThermoEnergy(meanRho, meanRhoU, meanE);
    const scalar meanRhoThermoEnergyFloor =
        calcRhoThermoEnergyFloor(meanRho);
    const scalar meanRhoThermoEnergyMargin =
        meanRhoThermoEnergy - meanRhoThermoEnergyFloor;

    theta2 =
        thetaCoeff
        (
            meanRhoThermoEnergyMargin,
            minRhoThermoEnergyMargin,
            min(scalar(0), meanRhoThermoEnergyMargin)
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
    thermoPtr_(nullptr),
    rhoMin_(dict.lookupOrDefault<scalar>("rhoMin", epsilon_)),
    rhoMax_(dict.lookupOrDefault<scalar>("rhoMax", GREAT)),
    magUMax_(dict.lookupOrDefault<scalar>("magUMax", GREAT)),
    Tmin_(dict.lookupOrDefault<scalar>("Tmin", epsilon_)),
    Tmax_(dict.lookupOrDefault<scalar>("Tmax", GREAT)),
    theta1_(1.0),
    theta2_(1.0),
    writeTheta_(dict.lookupOrDefault<bool>("writeTheta", false)),
    theta1FieldPtr_(nullptr),
    theta2FieldPtr_(nullptr)
{
    initializeConservativeFields();
    initializeThermo();

    if (epsilon_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'tolerance' must be positive, but got "
            << epsilon_ << '.'
            << exit(FatalIOError);
    }

    if (rhoMin_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'rhoMin' must be positive, but got "
            << rhoMin_ << '.'
            << exit(FatalIOError);
    }

    if (rhoMin_ > rhoMax_)
    {
        FatalIOErrorInFunction(dict)
            << "Entries 'rhoMin' and 'rhoMax' are inconsistent: rhoMin="
            << rhoMin_ << ", rhoMax=" << rhoMax_ << '.'
            << exit(FatalIOError);
    }

    if (magUMax_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'magUMax' must be non-negative, but got "
            << magUMax_ << '.'
            << exit(FatalIOError);
    }

    if (Tmin_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'Tmin' must be positive, but got "
            << Tmin_ << '.'
            << exit(FatalIOError);
    }

    if (Tmin_ > Tmax_)
    {
        FatalIOErrorInFunction(dict)
            << "Entries 'Tmin' and 'Tmax' are inconsistent: Tmin="
            << Tmin_ << ", Tmax=" << Tmax_ << '.'
            << exit(FatalIOError);
    }

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
    rhoMin_ = dict.lookupOrDefault<scalar>("rhoMin", epsilon_);
    rhoMax_ = dict.lookupOrDefault<scalar>("rhoMax", GREAT);
    magUMax_ = dict.lookupOrDefault<scalar>("magUMax", GREAT);
    Tmin_ = dict.lookupOrDefault<scalar>("Tmin", epsilon_);
    Tmax_ = dict.lookupOrDefault<scalar>("Tmax", GREAT);
    theta1_ = 1.0;
    theta2_ = 1.0;
    writeTheta_ = dict.lookupOrDefault<bool>("writeTheta", false);
    theta1FieldPtr_.clear();
    theta2FieldPtr_.clear();
    initializeConservativeFields();
    initializeThermo();

    if (epsilon_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'tolerance' must be positive, but got "
            << epsilon_ << '.'
            << exit(FatalIOError);
    }

    if (rhoMin_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'rhoMin' must be positive, but got "
            << rhoMin_ << '.'
            << exit(FatalIOError);
    }

    if (rhoMin_ > rhoMax_)
    {
        FatalIOErrorInFunction(dict)
            << "Entries 'rhoMin' and 'rhoMax' are inconsistent: rhoMin="
            << rhoMin_ << ", rhoMax=" << rhoMax_ << '.'
            << exit(FatalIOError);
    }

    if (magUMax_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'magUMax' must be non-negative, but got "
            << magUMax_ << '.'
            << exit(FatalIOError);
    }

    if (Tmin_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'Tmin' must be positive, but got "
            << Tmin_ << '.'
            << exit(FatalIOError);
    }

    if (Tmin_ > Tmax_)
    {
        FatalIOErrorInFunction(dict)
            << "Entries 'Tmin' and 'Tmax' are inconsistent: Tmin="
            << Tmin_ << ", Tmax=" << Tmax_ << '.'
            << exit(FatalIOError);
    }

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
        const bool meanLimited = limitCellMeanValue(cellID);

        if (meanLimited)
        {
            theta1_ = 0.0;
            theta2_ = 0.0;
        }
        else if
        (
            mesh_.pOrder() > 0
         && densityFieldPtr_->hasDof()
         && densityFieldPtr_->dof()[cellID].nDof() > 1
        )
        {
            computeThetaCoeffs(cellID, theta1_, theta2_);
        }
        else
        {
            // For p=0 only the mean remains, so skip the positivity scaling.
            theta1_ = 1.0;
            theta2_ = 1.0;
        }

        setThetaFields(cellID, theta1_, theta2_);

        const bool thetaLimited =
            meanLimited || theta1_ < 1.0 || theta2_ < 1.0;

        if (thetaLimited)
        {
            forAll(limitedFields_, fieldI)
            {
                limitCell(limitedFields_[fieldI], cellID);
            }
        }

        if (!thetaLimited)
        {
            continue;
        }

        ++nLimitedCells_;
        cacheLimitedCell(cellID);
    }

    postCorrect();
}

} // End namespace Foam

// ************************************************************************* //
