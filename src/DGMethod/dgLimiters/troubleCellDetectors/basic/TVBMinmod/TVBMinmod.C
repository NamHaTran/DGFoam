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

#include "TVBMinmod.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(TVBMinmod, 0);
addToRunTimeSelectionTable(troubleCellDetector, TVBMinmod, dictionary);

static scalar scalarMinmod(const List<scalar>& args);
static scalar tvbMinmodValue(const List<scalar>& args, const scalar condition);
static scalar vectorComponentValue(const vector& v, const direction cmpt);
static scalar cellMean(const dgField<scalar>& field, const label cellID);
static scalar cellMean
(
    const dgField<vector>& field,
    const label cellID,
    const direction cmpt
);
static scalar weightedFaceMean
(
    const faceGaussField<scalar>& faceField,
    const label localFaceI,
    const bool usePlus
);
static scalar weightedFaceMean
(
    const faceGaussField<vector>& faceField,
    const label localFaceI,
    const direction cmpt,
    const bool usePlus
);
static bool detectScalarField
(
    const dgGeomCell& cell,
    const label cellID,
    const scalar M,
    const bool includeBoundaryFaces,
    const dgField<scalar>& field,
    const GaussField<scalar>& gaussField
);
static bool detectVectorComponent
(
    const dgGeomCell& cell,
    const label cellID,
    const scalar M,
    const bool includeBoundaryFaces,
    const dgField<vector>& field,
    const GaussField<vector>& gaussField,
    const direction cmpt
);


TVBMinmod::TVBMinmod
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    troubleCellDetector(dict, mesh),
    M_(dict.lookupOrDefault<scalar>("M", 1.0)),
    includeBoundaryFaces_(dict.lookupOrDefault<bool>("includeBoundaryFaces", true))
{}


bool TVBMinmod::detect(const label cellID) const
{
    // Step 1: gather geometry and face/cell data for the target cell.
    const dgGeomCell& cell = *mesh_.cells()[cellID];

    // Step 2: loop through the configured checkFields list in order.
    forAll(checkFields(), fieldI)
    {
        const checkedField& checkField = checkFields()[fieldI];

        if (!checkField.isVector)
        {
            if
            (
                detectScalarField
                (
                    cell,
                    cellID,
                    M_,
                    includeBoundaryFaces_,
                    *checkField.scalarFieldPtr,
                    checkField.scalarFieldPtr->gaussFields()[cellID]
                )
            )
            {
                setLimitingIndicator(cellID, true);
                return true;
            }

            continue;
        }

        // Step 3: vector fields are checked component-by-component.
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            if
            (
                detectVectorComponent
                (
                    cell,
                    cellID,
                    M_,
                    includeBoundaryFaces_,
                    *checkField.vectorFieldPtr,
                    checkField.vectorFieldPtr->gaussFields()[cellID],
                    cmpt
                )
            )
            {
                setLimitingIndicator(cellID, true);
                return true;
            }
        }
    }

    setLimitingIndicator(cellID, false);
    return false;
}

/**
 * \brief Classical minmod operator for a list of scalar arguments.
 *
 * \param args Candidate values to be compared
 * \return The common-sign minimum magnitude, or zero if signs disagree
 */
static scalar scalarMinmod(const List<scalar>& args)
{
    if (args.empty())
    {
        return Zero;
    }

    const scalar sign0 = sign(args[0]);
    scalar minAbs = mag(args[0]);

    forAll(args, i)
    {
        if (sign(args[i]) != sign0)
        {
            return Zero;
        }

        minAbs = min(minAbs, mag(args[i]));
    }

    return sign0*minAbs;
}


/**
 * \brief TVB-relaxed minmod operator.
 *
 * \param args Input list whose first entry is the owner-face jump
 * \param condition TVB relaxation threshold \f$Mh^2\f$
 * \return The relaxed jump value
 */
static scalar tvbMinmodValue(const List<scalar>& args, const scalar condition)
{
    if (args.empty())
    {
        return Zero;
    }

    if (mag(args[0]) <= condition)
    {
        return args[0];
    }

    return scalarMinmod(args);
}


/**
 * \brief Extract one Cartesian component from a vector.
 *
 * \param v Input vector
 * \param cmpt Requested Cartesian component
 * \return Scalar value of the selected component
 */
static scalar vectorComponentValue(const vector& v, const direction cmpt)
{
    switch (cmpt)
    {
        case vector::X:
            return v.x();

        case vector::Y:
            return v.y();

        default:
            return v.z();
    }
}


/**
 * \brief Return the modal cell-average of a scalar DG field.
 *
 * \param field Scalar DG field
 * \param cellID Target cell index
 * \return Cell-average value stored in the zeroth modal coefficient
 */
static scalar cellMean(const dgField<scalar>& field, const label cellID)
{
    const List<scalar>& cellDof = field.dof()[cellID].dof();

    if (cellDof.empty())
    {
        FatalErrorInFunction
            << "Scalar DG field '" << field.name() << "' has no DoF in cell "
            << cellID << '.'
            << abort(FatalError);
    }

    return cellDof[0];
}


/**
 * \brief Return the modal cell-average of one vector component.
 *
 * \param field Vector DG field
 * \param cellID Target cell index
 * \param cmpt Requested Cartesian component
 * \return Mean value of the selected component in the target cell
 */
static scalar cellMean
(
    const dgField<vector>& field,
    const label cellID,
    const direction cmpt
)
{
    const List<vector>& cellDof = field.dof()[cellID].dof();

    if (cellDof.empty())
    {
        FatalErrorInFunction
            << "Vector DG field '" << field.name() << "' has no DoF in cell "
            << cellID << '.'
            << abort(FatalError);
    }

    return vectorComponentValue(cellDof[0], cmpt);
}


/**
 * \brief Compute a quadrature-weighted mean on one scalar face trace.
 *
 * \param faceField Face Gauss field
 * \param localFaceI Local face index inside the current cell
 * \param usePlus Select exterior trace when true, owner trace otherwise
 * \return Weighted face-average value
 */
static scalar weightedFaceMean
(
    const faceGaussField<scalar>& faceField,
    const label localFaceI,
    const bool usePlus
)
{
    const dgGeomFace& face = *faceField.faces()[localFaceI];
    const List<scalar>& weights = face.weights();
    const List<scalar>& J2D = face.J2D();

    scalar sumVal = Zero;
    scalar sumArea = Zero;

    forAll(weights, gpI)
    {
        const scalar ds = weights[gpI]*J2D[gpI];
        const scalar value =
            usePlus
          ? faceField.plusValueOnFace(localFaceI, gpI)
          : faceField.minusValueOnFace(localFaceI, gpI);

        sumVal += ds*value;
        sumArea += ds;
    }

    return sumVal/max(sumArea, VSMALL);
}


/**
 * \brief Compute a quadrature-weighted mean on one vector face trace.
 *
 * \param faceField Face Gauss field
 * \param localFaceI Local face index inside the current cell
 * \param cmpt Requested Cartesian component
 * \param usePlus Select exterior trace when true, owner trace otherwise
 * \return Weighted face-average of the selected vector component
 */
static scalar weightedFaceMean
(
    const faceGaussField<vector>& faceField,
    const label localFaceI,
    const direction cmpt,
    const bool usePlus
)
{
    const dgGeomFace& face = *faceField.faces()[localFaceI];
    const List<scalar>& weights = face.weights();
    const List<scalar>& J2D = face.J2D();

    scalar sumVal = Zero;
    scalar sumArea = Zero;

    forAll(weights, gpI)
    {
        const scalar ds = weights[gpI]*J2D[gpI];
        const vector value =
            usePlus
          ? faceField.plusValueOnFace(localFaceI, gpI)
          : faceField.minusValueOnFace(localFaceI, gpI);

        sumVal += ds*vectorComponentValue(value, cmpt);
        sumArea += ds;
    }

    return sumVal/max(sumArea, VSMALL);
}


/**
 * \brief Apply the face-by-face TVB-minmod check to one scalar DG field.
 *
 * \param cell Geometric cell providing size and connectivity
 * \param cellID Target cell index
 * \param M TVB constant in the relaxation threshold
 * \param includeBoundaryFaces Switch to evaluate boundary faces as well
 * \param field Scalar DG field to be checked
 * \param gaussField Cell-local Gauss field for the scalar variable
 * \return True if at least one face violates the TVB-minmod check
 */
static bool detectScalarField
(
    const dgGeomCell& cell,
    const label cellID,
    const scalar M,
    const bool includeBoundaryFaces,
    const dgField<scalar>& field,
    const GaussField<scalar>& gaussField
)
{
    const scalar meanP = cellMean(field, cellID);
    const scalar h = max(cell.cellSize(), VSMALL);
    const scalar condition = M*sqr(h);
    const labelList& neighbours = cell.neighborCells();
    const faceGaussField<scalar>& faceField = gaussField.faceField();

    List<scalar> args(2, Zero);

    forAll(neighbours, localFaceI)
    {
        const label neighID = neighbours[localFaceI];

        if (neighID < 0 && !includeBoundaryFaces)
        {
            continue;
        }

        args[0] = weightedFaceMean(faceField, localFaceI, false) - meanP;

        if (neighID >= 0)
        {
            args[1] = cellMean(field, neighID) - meanP;
        }
        else
        {
            args[1] = weightedFaceMean(faceField, localFaceI, true) - meanP;
        }

        if (mag(tvbMinmodValue(args, condition) - args[0]) > 1e-10)
        {
            return true;
        }
    }

    return false;
}


/**
 * \brief Apply the face-by-face TVB-minmod check to one vector component.
 *
 * \param cell Geometric cell providing size and connectivity
 * \param cellID Target cell index
 * \param M TVB constant in the relaxation threshold
 * \param includeBoundaryFaces Switch to evaluate boundary faces as well
 * \param field Vector DG field to be checked
 * \param gaussField Cell-local Gauss field for the vector variable
 * \param cmpt Requested Cartesian component
 * \return True if at least one face violates the TVB-minmod check
 */
static bool detectVectorComponent
(
    const dgGeomCell& cell,
    const label cellID,
    const scalar M,
    const bool includeBoundaryFaces,
    const dgField<vector>& field,
    const GaussField<vector>& gaussField,
    const direction cmpt
)
{
    const scalar meanP = cellMean(field, cellID, cmpt);
    const scalar h = max(cell.cellSize(), VSMALL);
    const scalar condition = M*sqr(h);
    const labelList& neighbours = cell.neighborCells();
    const faceGaussField<vector>& faceField = gaussField.faceField();

    List<scalar> args(2, Zero);

    forAll(neighbours, localFaceI)
    {
        const label neighID = neighbours[localFaceI];

        if (neighID < 0 && !includeBoundaryFaces)
        {
            continue;
        }

        args[0] = weightedFaceMean(faceField, localFaceI, cmpt, false) - meanP;

        if (neighID >= 0)
        {
            args[1] = cellMean(field, neighID, cmpt) - meanP;
        }
        else
        {
            args[1] = weightedFaceMean(faceField, localFaceI, cmpt, true) - meanP;
        }

        if (mag(tvbMinmodValue(args, condition) - args[0]) > 1e-10)
        {
            return true;
        }
    }

    return false;
}
} // End namespace Foam

// ************************************************************************* //
