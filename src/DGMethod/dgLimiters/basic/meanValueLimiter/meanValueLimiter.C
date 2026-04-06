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

#include "meanValueLimiter.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "DynamicList.H"
#include "error.H"

namespace Foam
{

defineTypeNameAndDebug(meanValueLimiter, 0);
addToRunTimeSelectionTable(dgLimiter, meanValueLimiter, dictionary);

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


dictionary normalizeMeanValueLimiterDict(const dictionary& inputDict)
{
    if (!inputDict.found("density"))
    {
        FatalIOErrorInFunction(inputDict)
            << "Missing required entry 'density' in "
            << "meanValue limiter dictionary."
            << exit(FatalIOError);
    }

    if (!inputDict.found("momentum"))
    {
        FatalIOErrorInFunction(inputDict)
            << "Missing required entry 'momentum' in "
            << "meanValue limiter dictionary."
            << exit(FatalIOError);
    }

    if (!inputDict.found("energy"))
    {
        FatalIOErrorInFunction(inputDict)
            << "Missing required entry 'energy' in "
            << "meanValue limiter dictionary."
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


void meanValueLimiter::initializeConservativeFields()
{
    densityFieldPtr_ = nullptr;
    momentumFieldPtr_ = nullptr;
    energyFieldPtr_ = nullptr;

    densityFieldName_ = dict_.get<word>("density");
    momentumFieldName_ = dict_.get<word>("momentum");
    energyFieldName_ = dict_.get<word>("energy");

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


void meanValueLimiter::initializeThermo()
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


scalar meanValueLimiter::safeDensityForDivision(const scalar rho) const
{
    return
        mag(rho) > epsilon_
      ? rho
      : (rho >= 0 ? epsilon_ : -epsilon_);
}


scalar meanValueLimiter::calcMeanTemperature
(
    const scalar rho,
    const vector& rhoU,
    const scalar E
) const
{
    const scalar rhoSafe = safeDensityForDivision(rho);
    const vector U = rhoU/rhoSafe;
    const scalar he = E/rhoSafe - 0.5*magSqr(U);

    return thermoPtr_->energyModel().calcTfromHe(he);
}


scalar meanValueLimiter::calcMeanTotalEnergy
(
    const scalar rho,
    const vector& U,
    const scalar T
) const
{
    const scalar he = thermoPtr_->energyModel().calcHe(T);
    scalar E = rho*(he + 0.5*magSqr(U));

    if (thermoPtr_->heIsEnthalpy())
    {
        E -= thermoPtr_->eos().calcPFromRhoT(rho, T);
    }

    return E;
}


void meanValueLimiter::zeroHighOrderModes
(
    dgField<scalar>& field,
    const label cellID
) const
{
    if (!field.hasDof())
    {
        return;
    }

    cellDof<scalar>& cellModes = field.dof()[cellID];

    for (label modeI = 1; modeI < cellModes.nDof(); ++modeI)
    {
        cellModes[modeI] = Zero;
    }
}


void meanValueLimiter::zeroHighOrderModes
(
    dgField<vector>& field,
    const label cellID
) const
{
    if (!field.hasDof())
    {
        return;
    }

    cellDof<vector>& cellModes = field.dof()[cellID];

    for (label modeI = 1; modeI < cellModes.nDof(); ++modeI)
    {
        cellModes[modeI] = Zero;
    }
}


bool meanValueLimiter::limitCellMeanValue(const label cellID)
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

    bool updateRho = false;
    bool updateRhoU = false;
    bool updateE = false;

    const scalar rho0NoClamp = rhoModes[0];
    const scalar rho0Clamped = min(max(rho0NoClamp, rhoMin_), rhoMax_);

    if (rho0Clamped != rho0NoClamp)
    {
        rhoModes[0] = rho0Clamped;
        zeroHighOrderModes(*densityFieldPtr_, cellID);
        updateRho = true;
    }

    const vector U0 = rhoUModes[0]/safeDensityForDivision(rho0NoClamp);
    const scalar magU0 = mag(U0);

    if (magU0 > magUMax_ && magU0 > SMALL)
    {
        const vector UClamped = U0*(magUMax_/magU0);
        rhoUModes[0] = rhoModes[0]*UClamped;
        zeroHighOrderModes(*momentumFieldPtr_, cellID);
        updateRhoU = true;
    }

    scalar T0 = calcMeanTemperature(rhoModes[0], rhoUModes[0], EModes[0]);
    const scalar TClamped = min(max(T0, Tmin_), Tmax_);

    if (TClamped != T0)
    {
        const vector UMean = rhoUModes[0]/safeDensityForDivision(rhoModes[0]);
        EModes[0] = calcMeanTotalEnergy(rhoModes[0], UMean, TClamped);
        zeroHighOrderModes(*energyFieldPtr_, cellID);
        updateE = true;
    }

    if (updateRho)
    {
        densityFieldPtr_->dof().updateCellDof(cellID);
    }

    if (updateRhoU)
    {
        momentumFieldPtr_->dof().updateCellDof(cellID);
    }

    if (updateE)
    {
        energyFieldPtr_->dof().updateCellDof(cellID);
    }

    return updateRho || updateRhoU || updateE;
}


void meanValueLimiter::postCorrect()
{
    dgLimiter::postCorrect();

    label nLimitedCells = nLimitedCells_;
    reduce(nLimitedCells, sumOp<label>());

    Info<< "Limiter " << type()
        << " applied on " << nLimitedCells
        << " cell(s)." << nl;
}


void meanValueLimiter::limitField
(
    dgField<scalar>& field,
    const label cellID
)
{
    (void)field;
    (void)cellID;
}


void meanValueLimiter::limitField
(
    dgField<vector>& field,
    const label cellID
)
{
    (void)field;
    (void)cellID;
}


meanValueLimiter::meanValueLimiter
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    dgLimiter(normalizeMeanValueLimiterDict(dict), mesh),
    densityFieldName_(word::null),
    momentumFieldName_(word::null),
    energyFieldName_(word::null),
    densityFieldPtr_(nullptr),
    momentumFieldPtr_(nullptr),
    energyFieldPtr_(nullptr),
    thermoPtr_(nullptr),
    epsilon_(dict.lookupOrDefault<scalar>("tolerance", 1.0e-12)),
    rhoMin_(dict.get<scalar>("rhoMin")),
    rhoMax_(dict.get<scalar>("rhoMax")),
    magUMax_(dict.get<scalar>("magUMax")),
    Tmin_(dict.get<scalar>("Tmin")),
    Tmax_(dict.get<scalar>("Tmax"))
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

    if (Tmin_ > Tmax_)
    {
        FatalIOErrorInFunction(dict)
            << "Entries 'Tmin' and 'Tmax' are inconsistent: Tmin="
            << Tmin_ << ", Tmax=" << Tmax_ << '.'
            << exit(FatalIOError);
    }
}


void meanValueLimiter::read(const dictionary& dict)
{
    dgLimiter::read(normalizeMeanValueLimiterDict(dict));
    epsilon_ = dict.lookupOrDefault<scalar>("tolerance", 1.0e-12);
    rhoMin_ = dict.get<scalar>("rhoMin");
    rhoMax_ = dict.get<scalar>("rhoMax");
    magUMax_ = dict.get<scalar>("magUMax");
    Tmin_ = dict.get<scalar>("Tmin");
    Tmax_ = dict.get<scalar>("Tmax");

    initializeConservativeFields();
    initializeThermo();

    if (epsilon_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'tolerance' must be positive, but got "
            << epsilon_ << '.'
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

    if (Tmin_ > Tmax_)
    {
        FatalIOErrorInFunction(dict)
            << "Entries 'Tmin' and 'Tmax' are inconsistent: Tmin="
            << Tmin_ << ", Tmax=" << Tmax_ << '.'
            << exit(FatalIOError);
    }
}


void meanValueLimiter::correct()
{
    preCorrect();
    nLimitedCells_ = 0;

    for (label cellID = 0; cellID < mesh_.nCells(); ++cellID)
    {
        if (!limitCellMeanValue(cellID))
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
