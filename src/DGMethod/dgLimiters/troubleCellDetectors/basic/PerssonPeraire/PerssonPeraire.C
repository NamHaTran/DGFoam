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

#include "PerssonPeraire.H"
#include <cmath>
#include "DynamicList.H"
#include "addToRunTimeSelectionTable.H"
#include "dgCellType.H"
#include "error.H"

namespace Foam
{

namespace
{

scalar defaultS0(const label pOrder)
{
    const scalar Se0 = 1.0/std::pow(scalar(pOrder + 1), 4.0);
    return std::log10(max(Se0, VSMALL));
}

dgCellType inferCellType(const dgGeomCell& cell)
{
    switch (cell.nPoints())
    {
        case 4:
            return dgCellType::TET;

        case 5:
            return dgCellType::PYRAMID;

        case 6:
            return dgCellType::PRISM;

        case 8:
            return dgCellType::HEX;

        default:
            FatalErrorInFunction
                << "Unsupported cell topology with " << cell.nPoints()
                << " points and " << cell.nFaces() << " faces."
                << abort(FatalError);
    }

    return dgCellType::INVALID;
}


List<label> pOrderModeIndices
(
    const dgCellType cellType,
    const label pOrder,
    label& nDof
)
{
    DynamicList<label> indices;
    label modeI = 0;

    switch (cellType)
    {
        case dgCellType::HEX:
        case dgCellType::PYRAMID:
        case dgCellType::TET:
        {
            for (label i = 0; i <= pOrder; ++i)
            {
                for (label j = 0; j <= pOrder - i; ++j)
                {
                    for (label k = 0; k <= pOrder - i - j; ++k)
                    {
                        if (i + j + k == pOrder)
                        {
                            indices.append(modeI);
                        }

                        ++modeI;
                    }
                }
            }

            break;
        }

        case dgCellType::PRISM:
        {
            for (label i = 0; i <= pOrder; ++i)
            {
                for (label j = 0; j <= pOrder; ++j)
                {
                    for (label k = 0; k <= pOrder - j; ++k)
                    {
                        if (i == pOrder || j + k == pOrder)
                        {
                            indices.append(modeI);
                        }

                        ++modeI;
                    }
                }
            }

            break;
        }

        default:
            FatalErrorInFunction
                << "Unsupported cell type " << cellType
                << " when collecting p-order modal indices."
                << abort(FatalError);
    }

    nDof = modeI;

    List<label> cached(indices.size());

    forAll(indices, i)
    {
        cached[i] = indices[i];
    }

    return cached;
}


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

} // End anonymous namespace

defineTypeNameAndDebug(PerssonPeraire, 0);
addToRunTimeSelectionTable(troubleCellDetector, PerssonPeraire, dictionary);


PerssonPeraire::PerssonPeraire
(
    const dictionary& dict,
    const dgGeomMesh& mesh
)
:
    troubleCellDetector(dict, mesh),
    tetNDof_(Zero),
    hexNDof_(Zero),
    prismNDof_(Zero),
    pyramidNDof_(Zero),
    tetPModeDof_(pOrderModeIndices(dgCellType::TET, mesh.pOrder(), tetNDof_)),
    hexPModeDof_(pOrderModeIndices(dgCellType::HEX, mesh.pOrder(), hexNDof_)),
    prismPModeDof_(pOrderModeIndices(dgCellType::PRISM, mesh.pOrder(), prismNDof_)),
    pyramidPModeDof_(pOrderModeIndices(dgCellType::PYRAMID, mesh.pOrder(), pyramidNDof_)),
    s0_(dict.lookupOrDefault<scalar>("s0", defaultS0(mesh.pOrder()))),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 1.0))
{}


const List<label>& PerssonPeraire::pModeDof(const dgCellType cellType) const
{
    switch (cellType)
    {
        case dgCellType::TET:
            return tetPModeDof_;

        case dgCellType::HEX:
            return hexPModeDof_;

        case dgCellType::PRISM:
            return prismPModeDof_;

        case dgCellType::PYRAMID:
            return pyramidPModeDof_;

        default:
            FatalErrorInFunction
                << "Unsupported cell type " << cellType
                << " when accessing cached p-order shell indices."
                << abort(FatalError);
    }

    return tetPModeDof_;
}


label PerssonPeraire::expectedNDof(const dgCellType cellType) const
{
    switch (cellType)
    {
        case dgCellType::TET:
            return tetNDof_;

        case dgCellType::HEX:
            return hexNDof_;

        case dgCellType::PRISM:
            return prismNDof_;

        case dgCellType::PYRAMID:
            return pyramidNDof_;

        default:
            FatalErrorInFunction
                << "Unsupported cell type " << cellType
                << " when accessing cached modal size."
                << abort(FatalError);
    }

    return tetNDof_;
}


bool PerssonPeraire::detect(const label cellID) const
{
    const dgGeomCell& cell = *mesh_.cells()[cellID];
    const dgCellType cellType = inferCellType(cell);
    const List<label>& pModeDof = this->pModeDof(cellType);
    const label expectedDof = expectedNDof(cellType);
    const List<List<scalar>>& basis = cell.basis();
    const List<scalar>& weights = cell.weights();
    const List<scalar>& J3D = cell.J3D();

    forAll(checkFields(), fieldI)
    {
        const checkedField& checkField = checkFields()[fieldI];

        if (!checkField.isVector)
        {
            const dgField<scalar>& field = *checkField.scalarFieldPtr;
            const GaussField<scalar>& fieldGF = field.gaussFields()[cellID];
            const cellGaussField<scalar>& fieldCell = fieldGF.cellField();
            const dofField<scalar>* fieldDofPtr = fieldGF.dofFieldPtr();

            if (!fieldDofPtr)
            {
                FatalErrorInFunction
                    << "PerssonPeraire detector requires a scalar dofField for '"
                    << field.name() << "'."
                    << abort(FatalError);
            }

            const List<scalar>& fieldDof = (*fieldDofPtr)[cellID].dof();

            if (fieldDof.size() <= 1)
            {
                continue;
            }

            if (fieldDof.size() != expectedDof)
            {
                FatalErrorInFunction
                    << "Modal indexing mismatch for cell type " << cellType
                    << ": cached " << expectedDof
                    << " modes but scalar field '" << field.name()
                    << "' stores " << fieldDof.size() << '.'
                    << abort(FatalError);
            }

            scalar numerator = Zero;
            scalar denominator = Zero;

            forAll(weights, gpI)
            {
                scalar highOrderVal = Zero;

                forAll(pModeDof, shellI)
                {
                    const label dofI = pModeDof[shellI];
                    highOrderVal += basis[gpI][dofI]*fieldDof[dofI];
                }

                const scalar w = weights[gpI]*J3D[gpI];
                const scalar fieldVal = fieldCell[gpI];

                numerator += w*sqr(highOrderVal);
                denominator += w*sqr(fieldVal);
            }

            if (denominator <= VSMALL)
            {
                continue;
            }

            const scalar sensor = numerator/denominator;
            const scalar logSensor = std::log10(max(sensor, VSMALL));

            if (logSensor >= (s0_ - kappa_))
            {
                setLimitingIndicator(cellID, true);
                return true;
            }

            continue;
        }

        const dgField<vector>& field = *checkField.vectorFieldPtr;
        const GaussField<vector>& fieldGF = field.gaussFields()[cellID];
        const cellGaussField<vector>& fieldCell = fieldGF.cellField();
        const dofField<vector>* fieldDofPtr = fieldGF.dofFieldPtr();

        if (!fieldDofPtr)
        {
            FatalErrorInFunction
                << "PerssonPeraire detector requires a vector dofField for '"
                << field.name() << "'."
                << abort(FatalError);
        }

        const List<vector>& fieldDof = (*fieldDofPtr)[cellID].dof();

        if (fieldDof.size() <= 1)
        {
            continue;
        }

        if (fieldDof.size() != expectedDof)
        {
            FatalErrorInFunction
                << "Modal indexing mismatch for cell type " << cellType
                << ": cached " << expectedDof
                << " modes but vector field '" << field.name()
                << "' stores " << fieldDof.size() << '.'
                << abort(FatalError);
        }

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            scalar numerator = Zero;
            scalar denominator = Zero;

            forAll(weights, gpI)
            {
                scalar highOrderVal = Zero;

                forAll(pModeDof, shellI)
                {
                    const label dofI = pModeDof[shellI];
                    highOrderVal +=
                        basis[gpI][dofI]
                       *vectorComponentValue(fieldDof[dofI], cmpt);
                }

                const scalar w = weights[gpI]*J3D[gpI];
                const scalar fieldVal =
                    vectorComponentValue(fieldCell[gpI], cmpt);

                numerator += w*sqr(highOrderVal);
                denominator += w*sqr(fieldVal);
            }

            if (denominator <= VSMALL)
            {
                continue;
            }

            const scalar sensor = numerator/denominator;
            const scalar logSensor = std::log10(max(sensor, VSMALL));

            if (logSensor >= (s0_ - kappa_))
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
