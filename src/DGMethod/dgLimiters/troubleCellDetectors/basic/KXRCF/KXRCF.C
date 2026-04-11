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

#include "KXRCF.H"
#include <cmath>
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace
{

scalar cellMaxAbs(const GaussField<scalar>& gf)
{
    const cellGaussField<scalar>& cellField = gf.cellField();
    scalar maxVal = Zero;

    for (label gpI = 0; gpI < cellField.size(); ++gpI)
    {
        maxVal = max(maxVal, mag(cellField[gpI]));
    }

    return maxVal;
}


scalar cellMaxAbs
(
    const GaussField<vector>& gf,
    const direction cmpt
)
{
    const cellGaussField<vector>& cellField = gf.cellField();
    scalar maxVal = Zero;

    for (label gpI = 0; gpI < cellField.size(); ++gpI)
    {
        maxVal = max(maxVal, mag(cellField[gpI].component(cmpt)));
    }

    return maxVal;
}


scalar cellMean
(
    const dgField<scalar>& field,
    const label cellID
)
{
    const dofField<scalar>& dof = field.dof();
    const List<scalar>& cellDof = dof[cellID].dof();

    if (cellDof.empty())
    {
        FatalErrorInFunction
            << "Scalar DG field '" << field.name() << "' has no DoF in cell "
            << cellID << '.'
            << abort(FatalError);
    }

    // In modal DG, the zeroth coefficient is the cell average.
    return cellDof[0];
}


scalar cellMean
(
    const dgField<vector>& field,
    const label cellID,
    const direction cmpt
)
{
    const dofField<vector>& dof = field.dof();
    const List<vector>& cellDof = dof[cellID].dof();

    if (cellDof.empty())
    {
        FatalErrorInFunction
            << "Vector DG field '" << field.name() << "' has no DoF in cell "
            << cellID << '.'
            << abort(FatalError);
    }

    return cellDof[0].component(cmpt);
}

} // End anonymous namespace

defineTypeNameAndDebug(KXRCF, 0);
addToRunTimeSelectionTable(troubleCellDetector, KXRCF, dictionary);


KXRCF::KXRCF
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    troubleCellDetector(dict, mesh),
    directionField_
    (
        mesh.getFvMesh().lookupObject<dgField<vector>>
        (
            lookupFieldName(dict, "directionField", "rhoU")
        )
    ),
    LPRMode_(dict.lookupOrDefault<bool>("LPRMode", false)),
    threshold_(dict.lookupOrDefault<scalar>("threshold", 1.0)),
    includeBoundaryFaces_
    (
        dict.lookupOrDefault<bool>("includeBoundaryFaces", true)
    )
{}


bool KXRCF::detect(const label cellID) const
{
    const dgGeomCell& cell = *mesh_.cells()[cellID];
    const GaussField<vector>& directionGF = directionField_.gaussFields()[cellID];
    const faceGaussField<vector>& directionFace = directionGF.faceField();

    forAll(checkFields(), fieldI)
    {
        const checkedField& checkField = checkFields()[fieldI];

        if (!checkField.isVector)
        {
            const dgField<scalar>& qField = *checkField.scalarFieldPtr;
            const GaussField<scalar>& qGF = qField.gaussFields()[cellID];
            const faceGaussField<scalar>& qFace = qGF.faceField();

            scalar inflowJumpIntegral = Zero;
            scalar inflowMeasure = Zero;
            scalar qNorm = cellMaxAbs(qGF);
            const scalar qMean = mag(cellMean(qField, cellID));

            for (label localFaceI = 0; localFaceI < qFace.nFaces(); ++localFaceI)
            {
                const dgGeomFace& face = *qFace.faces()[localFaceI];

                if (!includeBoundaryFaces_ && face.isBoundary())
                {
                    continue;
                }

                const List<scalar>& weights = face.weights();
                const List<scalar>& J2D = face.J2D();

                forAll(weights, gpI)
                {
                    const vector directionVal =
                        directionFace.minusValueOnFace(localFaceI, gpI);
                    const scalar beta =
                        directionVal & directionFace.normals()[localFaceI];

                    const scalar qMinus = qFace.minusValueOnFace(localFaceI, gpI);
                    qNorm = max(qNorm, mag(qMinus));

                    if (beta >= 0)
                    {
                        continue;
                    }

                    const scalar ds = weights[gpI]*J2D[gpI];
                    const scalar qPlus = qFace.plusValueOnFace(localFaceI, gpI);

                    inflowJumpIntegral += (qMinus - qPlus)*ds;
                    inflowMeasure += ds;
                }
            }

            if (inflowMeasure <= VSMALL || qNorm <= VSMALL)
            {
                continue;
            }

            scalar sensor = Zero;

            if (LPRMode_)
            {
                const scalar hScale = max(cell.cellSize(), VSMALL);
                sensor =
                    mag(inflowJumpIntegral)
                   /(hScale*inflowMeasure*max(qMean, VSMALL));
            }
            else
            {
                const scalar hScale =
                    std::pow(max(cell.cellSize(), VSMALL), 0.5*scalar(mesh_.pOrder() + 1));

                sensor = mag(inflowJumpIntegral)/(hScale*inflowMeasure*qNorm);
            }

            if (sensor > threshold_)
            {
                setLimitingIndicator(cellID, true);
                return true;
            }

            continue;
        }

        const dgField<vector>& qField = *checkField.vectorFieldPtr;
        const GaussField<vector>& qGF = qField.gaussFields()[cellID];
        const faceGaussField<vector>& qFace = qGF.faceField();

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            scalar inflowJumpIntegral = Zero;
            scalar inflowMeasure = Zero;
            scalar qNorm = cellMaxAbs(qGF, cmpt);
            const scalar qMean = mag(cellMean(qField, cellID, cmpt));

            for (label localFaceI = 0; localFaceI < qFace.nFaces(); ++localFaceI)
            {
                const dgGeomFace& face = *qFace.faces()[localFaceI];

                if (!includeBoundaryFaces_ && face.isBoundary())
                {
                    continue;
                }

                const List<scalar>& weights = face.weights();
                const List<scalar>& J2D = face.J2D();

                forAll(weights, gpI)
                {
                    const vector directionVal =
                        directionFace.minusValueOnFace(localFaceI, gpI);
                    const scalar beta =
                        directionVal & directionFace.normals()[localFaceI];

                    const scalar qMinus =
                        qFace.minusValueOnFace(localFaceI, gpI).component(cmpt);
                    qNorm = max(qNorm, mag(qMinus));

                    if (beta >= 0)
                    {
                        continue;
                    }

                    const scalar ds = weights[gpI]*J2D[gpI];
                    const scalar qPlus =
                        qFace.plusValueOnFace(localFaceI, gpI).component(cmpt);

                    inflowJumpIntegral += (qMinus - qPlus)*ds;
                    inflowMeasure += ds;
                }
            }

            if (inflowMeasure <= VSMALL || qNorm <= VSMALL)
            {
                continue;
            }

            scalar sensor = Zero;

            if (LPRMode_)
            {
                const scalar hScale = max(cell.cellSize(), VSMALL);
                sensor =
                    mag(inflowJumpIntegral)
                   /(hScale*inflowMeasure*max(qMean, VSMALL));
            }
            else
            {
                const scalar hScale =
                    std::pow(max(cell.cellSize(), VSMALL), 0.5*scalar(mesh_.pOrder() + 1));

                sensor = mag(inflowJumpIntegral)/(hScale*inflowMeasure*qNorm);
            }

            if (sensor > threshold_)
            {
                setLimitingIndicator(cellID, true);
                return true;
            }
        }
    }

    setLimitingIndicator(cellID, false);
    return false;
}

} // End namespace Foam

// ************************************************************************* //
