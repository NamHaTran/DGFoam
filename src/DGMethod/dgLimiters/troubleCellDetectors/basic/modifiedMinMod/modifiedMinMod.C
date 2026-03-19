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

#include "modifiedMinMod.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace
{

/**
 * \brief Classical minmod function for a list of scalar arguments.
 *
 * \param args Input arguments to be compared
 * \return The minmod-limited value
 */
scalar scalarMinmod(const List<scalar>& args)
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
 * \brief TVB-style modified minmod used by the legacy detector.
 *
 * \param args Input arguments whose first entry is the candidate slope
 * \param condition Relaxation threshold, typically \f$M h^2 + \epsilon\f$
 * \return The relaxed minmod value
 */
scalar modifiedMinmodValue(const List<scalar>& args, const scalar condition)
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
 * \brief Extract one Cartesian component from a vector value.
 *
 * \param v Input vector
 * \param cmpt Requested component
 * \return Scalar value of the selected component
 */
scalar vectorComponentValue(const vector& v, const direction cmpt)
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
 * \brief Return the mean value of a scalar DG field in one cell.
 *
 * \param field Scalar DG field
 * \param cellID Target cell index
 * \return Mean cell value
 */
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
            << cellID << '.' << abort(FatalError);
    }

    // In modal DG, the cell average is the zeroth modal coefficient.
    return cellDof[0];
}


/**
 * \brief Return the mean value of one component of a vector DG field.
 *
 * \param field Vector DG field
 * \param cellID Target cell index
 * \param cmpt Requested Cartesian component
 * \return Mean component value in the cell
 */
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
            << cellID << '.' << abort(FatalError);
    }

    // In modal DG, the cell average is stored in the zeroth modal coefficient.
    return vectorComponentValue(cellDof[0], cmpt);
}


/**
 * \brief Approximate the owner-side face mean from face Gauss values.
 *
 * \param gf Cell-local scalar Gauss field
 * \param localFaceI Local face index inside the current cell
 * \return Arithmetic mean over owner-side face Gauss values
 */
scalar faceMean
(
    const GaussField<scalar>& gf,
    const label localFaceI
)
{
    const faceGaussField<scalar>& faceField = gf.faceField();
    scalar sumVal = Zero;

    for (label gpI = 0; gpI < faceField.nGaussPerFace(); ++gpI)
    {
        sumVal += faceField.minusValueOnFace(localFaceI, gpI);
    }

    return sumVal/max(label(1), faceField.nGaussPerFace());
}


/**
 * \brief Approximate the owner-side face mean of one vector component.
 *
 * \param gf Cell-local vector Gauss field
 * \param localFaceI Local face index inside the current cell
 * \param cmpt Requested Cartesian component
 * \return Arithmetic mean over owner-side face Gauss values
 */
scalar faceMean
(
    const GaussField<vector>& gf,
    const label localFaceI,
    const direction cmpt
)
{
    const faceGaussField<vector>& faceField = gf.faceField();
    scalar sumVal = Zero;

    for (label gpI = 0; gpI < faceField.nGaussPerFace(); ++gpI)
    {
        sumVal += vectorComponentValue
        (
            faceField.minusValueOnFace(localFaceI, gpI),
            cmpt
        );
    }

    return sumVal/max(label(1), faceField.nGaussPerFace());
}


/**
 * \brief Compute the legacy smoothness scale \f$M\f$ from neighboring means.
 *
 * \param cell Geometric cell providing neighbor connectivity
 * \param cellID Target cell index
 * \param meanFunctor Callable returning the mean value of a cell
 * \return Absolute value of the accumulated neighbor difference
 */
template<class MeanFunctor>
scalar calcM
(
    const dgGeomCell& cell,
    const label cellID,
    const MeanFunctor& meanFunctor
)
{
    const scalar meanP = meanFunctor(cellID);
    scalar M = Zero;

    const labelList& neighbours = cell.neighborCells();

    forAll(neighbours, localFaceI)
    {
        const label neighID = neighbours[localFaceI];
        const scalar neighMean =
            neighID >= 0 ? meanFunctor(neighID) : meanP;

        M += neighMean - meanP;
    }

    return mag(M);
}


/**
 * \brief Apply the modified-minmod detector to one scalar-like quantity.
 *
 * The detector is run face by face. For each internal face, the owner-side
 * face mean is compared against the set of neighbor-mean jumps using the
 * modified-minmod relaxation with threshold \f$M h^2 + 0.1\f$.
 *
 * \param cell Geometric cell providing neighbor connectivity and size
 * \param cellID Target cell index
 * \param meanFunctor Callable returning a cell mean
 * \param faceMeanFunctor Callable returning an owner-side face mean
 * \return True if the quantity is detected as troubled
 */
template<class MeanFunctor, class FaceMeanFunctor>
bool detectScalarLike
(
    const dgGeomCell& cell,
    const label cellID,
    const MeanFunctor& meanFunctor,
    const FaceMeanFunctor& faceMeanFunctor
)
{
    const scalar meanP = meanFunctor(cellID);
    const scalar M = calcM(cell, cellID, meanFunctor);
    const scalar condition = M*sqr(cell.cellSize()) + 0.1;
    const labelList& neighbours = cell.neighborCells();

    DynamicList<scalar> neighbourDeltas(neighbours.size());

    forAll(neighbours, localFaceI)
    {
        const label neighID = neighbours[localFaceI];

        if (neighID >= 0)
        {
            neighbourDeltas.append(meanFunctor(neighID) - meanP);
        }
    }

    if (neighbourDeltas.empty())
    {
        return false;
    }

    List<scalar> args(neighbourDeltas.size() + 1, Zero);

    forAll(neighbours, localFaceI)
    {
        if (neighbours[localFaceI] < 0)
        {
            continue;
        }

        args[0] = faceMeanFunctor(localFaceI) - meanP;

        forAll(neighbourDeltas, i)
        {
            args[i + 1] = neighbourDeltas[i];
        }

        if (mag(modifiedMinmodValue(args, condition) - args[0]) > 1e-10)
        {
            return true;
        }
    }

    return false;
}

} // End anonymous namespace


defineTypeNameAndDebug(modifiedMinMod, 0);
addToRunTimeSelectionTable(troubleCellDetector, modifiedMinMod, dictionary);


modifiedMinMod::modifiedMinMod
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    troubleCellDetector(dict, mesh)
{}


bool modifiedMinMod::detect(const label cellID) const
{
    const dgGeomCell& cell = *mesh_.cells()[cellID];

    const GaussField<scalar>& rhoGF = rhoGauss(cellID);
    const GaussField<vector>& rhoUGF = rhoUGauss(cellID);
    const GaussField<scalar>& EGF = EGauss(cellID);

    const auto rhoMean =
        [&](const label cI) -> scalar
        {
            return cellMean(rho_, cI);
        };

    const auto rhoFaceMean =
        [&](const label localFaceI) -> scalar
        {
            return faceMean(rhoGF, localFaceI);
        };

    // Mark the cell as troubled as soon as one conservative component
    // violates the modified-minmod check.
    if (detectScalarLike(cell, cellID, rhoMean, rhoFaceMean))
    {
        return true;
    }

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const auto rhoUMean =
            [&](const label cI) -> scalar
            {
                return cellMean(rhoU_, cI, cmpt);
            };

        const auto rhoUFaceMean =
            [&](const label localFaceI) -> scalar
            {
                return faceMean(rhoUGF, localFaceI, cmpt);
            };

        if (detectScalarLike(cell, cellID, rhoUMean, rhoUFaceMean))
        {
            return true;
        }
    }

    const auto EMean =
        [&](const label cI) -> scalar
        {
            return cellMean(E_, cI);
        };

    const auto EFaceMean =
        [&](const label localFaceI) -> scalar
        {
            return faceMean(EGF, localFaceI);
        };

    return detectScalarLike(cell, cellID, EMean, EFaceMean);
}

} // End namespace Foam

// ************************************************************************* //
