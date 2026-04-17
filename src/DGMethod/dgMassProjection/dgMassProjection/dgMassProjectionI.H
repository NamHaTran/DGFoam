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

#ifndef Foam_dgMassProjectionI_H
#define Foam_dgMassProjectionI_H

#include "dgMassProjection.H"
#include "PstreamReduceOps.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
inline void dgMassProjection::checkType() const
{
    static_assert
    (
        std::is_same_v<Type, scalar>
     || std::is_same_v<Type, vector>
     || std::is_same_v<Type, tensor>,
        "dgMassProjection: unsupported field type. "
        "Allowed: scalar, vector, tensor"
    );
}


template<class Type>
inline auto& dgMassProjection::names()
{
    checkType<Type>();

    if constexpr (std::is_same_v<Type, scalar>)
    {
        return scalarNames_;
    }
    else if constexpr (std::is_same_v<Type, vector>)
    {
        return vectorNames_;
    }
    else if constexpr (std::is_same_v<Type, tensor>)
    {
        return tensorNames_;
    }
    else
    {
        static_assert(alwaysFalse<Type>::value, "Unsupported field type");
    }
}


template<class Type>
inline const auto& dgMassProjection::names() const
{
    checkType<Type>();

    if constexpr (std::is_same_v<Type, scalar>)
    {
        return scalarNames_;
    }
    else if constexpr (std::is_same_v<Type, vector>)
    {
        return vectorNames_;
    }
    else if constexpr (std::is_same_v<Type, tensor>)
    {
        return tensorNames_;
    }
    else
    {
        static_assert(alwaysFalse<Type>::value, "Unsupported field type");
    }
}


template<class Type>
inline auto& dgMassProjection::fields()
{
    checkType<Type>();

    if constexpr (std::is_same_v<Type, scalar>)
    {
        return scalarFields_;
    }
    else if constexpr (std::is_same_v<Type, vector>)
    {
        return vectorFields_;
    }
    else if constexpr (std::is_same_v<Type, tensor>)
    {
        return tensorFields_;
    }
    else
    {
        static_assert(alwaysFalse<Type>::value, "Unsupported field type");
    }
}


template<class Type>
inline const auto& dgMassProjection::fields() const
{
    checkType<Type>();

    if constexpr (std::is_same_v<Type, scalar>)
    {
        return scalarFields_;
    }
    else if constexpr (std::is_same_v<Type, vector>)
    {
        return vectorFields_;
    }
    else if constexpr (std::is_same_v<Type, tensor>)
    {
        return tensorFields_;
    }
    else
    {
        static_assert(alwaysFalse<Type>::value, "Unsupported field type");
    }
}


template<class Type>
inline auto& dgMassProjection::residuals()
{
    checkType<Type>();

    if constexpr (std::is_same_v<Type, scalar>)
    {
        return scalarResiduals_;
    }
    else if constexpr (std::is_same_v<Type, vector>)
    {
        return vectorResiduals_;
    }
    else if constexpr (std::is_same_v<Type, tensor>)
    {
        return tensorResiduals_;
    }
    else
    {
        static_assert(alwaysFalse<Type>::value, "Unsupported field type");
    }
}


template<class Type>
inline const auto& dgMassProjection::residuals() const
{
    checkType<Type>();

    if constexpr (std::is_same_v<Type, scalar>)
    {
        return scalarResiduals_;
    }
    else if constexpr (std::is_same_v<Type, vector>)
    {
        return vectorResiduals_;
    }
    else if constexpr (std::is_same_v<Type, tensor>)
    {
        return tensorResiduals_;
    }
    else
    {
        static_assert(alwaysFalse<Type>::value, "Unsupported field type");
    }
}


template<class Type>
inline List<List<Type>>& dgMassProjection::residualBuffer
(
    const word& fieldName
)
{
    List<word>& fieldNames = names<Type>();
    const label fieldI = findNameIndex(fieldNames, fieldName);

    if (fieldI < 0)
    {
        FatalErrorInFunction
            << "Field " << fieldName
            << " is not registered in dgMassProjection as "
            << pTraits<Type>::typeName << '.'
            << abort(FatalError);
    }

    return residuals<Type>()[fieldI];
}


template<class Type>
inline const List<List<Type>>& dgMassProjection::residualBuffer
(
    const word& fieldName
) const
{
    const List<word>& fieldNames = names<Type>();
    const label fieldI = findNameIndex(fieldNames, fieldName);

    if (fieldI < 0)
    {
        FatalErrorInFunction
            << "Field " << fieldName
            << " is not registered in dgMassProjection as "
            << pTraits<Type>::typeName << '.'
            << abort(FatalError);
    }

    return residuals<Type>()[fieldI];
}


template<class Type>
void dgMassProjection::allocateResidual(List<List<Type>>& R) const
{
    R.setSize(mesh_.nCells());

    forAll(R, cellI)
    {
        const label nDof = mesh_.cells()[cellI]->nDof();
        R[cellI].setSize(nDof, pTraits<Type>::zero);
    }
}


template<class Type>
void dgMassProjection::clearResidualCell
(
    List<List<Type>>& R,
    const label cellID
) const
{
    validateCellID(cellID);

    forAll(R[cellID], dofI)
    {
        R[cellID][dofI] = pTraits<Type>::zero;
    }
}


template<class Type>
void dgMassProjection::clearResidual(List<List<Type>>& R) const
{
    forAll(R, cellI)
    {
        clearResidualCell(R, cellI);
    }
}


template<class Type>
void dgMassProjection::solveFieldCell
(
    dgField<Type>& field,
    const List<List<Type>>& R,
    const label cellID
) const
{
    validateCellID(cellID);

    if (!field.hasDof())
    {
        FatalErrorInFunction
            << "Cannot solve auxiliary field " << field.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    if (cellID >= R.size())
    {
        FatalErrorInFunction
            << "Residual storage for field " << field.name()
            << " has only " << R.size() << " cells; requested cell "
            << cellID << '.'
            << abort(FatalError);
    }

    List<Type>& cellDof = field.dof()[cellID].dof();
    const List<Type>& cellResidual = R[cellID];
    const List<scalar> massDiag(mesh_.cells()[cellID]->massMatrixDiag());

    if
    (
        cellDof.size() != cellResidual.size()
     || cellDof.size() != massDiag.size()
    )
    {
        FatalErrorInFunction
            << "Size mismatch while solving auxiliary field " << field.name()
            << " in cell " << cellID << ": dof=" << cellDof.size()
            << ", residual=" << cellResidual.size()
            << ", massDiag=" << massDiag.size() << '.'
            << abort(FatalError);
    }

    forAll(cellDof, dofI)
    {
        if (mag(massDiag[dofI]) <= VSMALL)
        {
            FatalErrorInFunction
                << "Near-zero mass-matrix diagonal for field " << field.name()
                << " in cell " << cellID << ", dof " << dofI
                << ": " << massDiag[dofI]
                << abort(FatalError);
        }

        cellDof[dofI] = cellResidual[dofI]/massDiag[dofI];
    }

    field.gaussFields()[cellID].interpolateFromDof();
    field.dof().updateCellDof(cellID);
}


template<class Type>
void dgMassProjection::solveField
(
    dgField<Type>& field,
    const List<List<Type>>& R
) const
{
    for (label cellID = 0; cellID < mesh_.nCells(); ++cellID)
    {
        solveFieldCell(field, R, cellID);
    }
}


template<class Type>
void dgMassProjection::solveType(const label cellID) const
{
    const auto& fieldPtrs = fields<Type>();
    const auto& fieldResiduals = residuals<Type>();

    forAll(fieldPtrs, fieldI)
    {
        solveFieldCell(*fieldPtrs[fieldI], fieldResiduals[fieldI], cellID);
    }
}


template<class Type>
void dgMassProjection::solveType() const
{
    const auto& fieldPtrs = fields<Type>();
    const auto& fieldResiduals = residuals<Type>();

    forAll(fieldPtrs, fieldI)
    {
        solveField(*fieldPtrs[fieldI], fieldResiduals[fieldI]);
    }
}


template<class Type>
void dgMassProjection::writeResidualType
(
    Ostream& os,
    const word& typeLabel
) const
{
    const auto& fieldNames = names<Type>();
    const auto& fieldResiduals = residuals<Type>();

    forAll(fieldNames, fieldI)
    {
        scalar maxResidual = Zero;

        forAll(fieldResiduals[fieldI], cellI)
        {
            const List<Type>& cellResidual = fieldResiduals[fieldI][cellI];

            forAll(cellResidual, dofI)
            {
                maxResidual = max(maxResidual, mag(cellResidual[dofI]));
            }
        }

        reduce(maxResidual, maxOp<scalar>());

        os  << "  " << typeLabel << ' ' << fieldNames[fieldI]
            << " : max(|R|) = " << maxResidual << nl;
    }
}


// * * * * * * * * * * * * * Public Template Functions * * * * * * * * * * * //

template<class Type>
void dgMassProjection::registerField(dgField<Type>& field)
{
    checkType<Type>();

    if (!field.hasDof())
    {
        FatalErrorInFunction
            << "Cannot register auxiliary field " << field.name()
            << " because it has no DOF storage."
            << abort(FatalError);
    }

    if (field.nCells() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Cannot register auxiliary field " << field.name()
            << " because it has " << field.nCells()
            << " cells but the DG mesh has " << mesh_.nCells() << '.'
            << abort(FatalError);
    }

    List<word>& fieldNames = names<Type>();

    if (findNameIndex(fieldNames, field.name()) >= 0)
    {
        FatalErrorInFunction
            << "Field " << field.name()
            << " is already registered in dgMassProjection."
            << abort(FatalError);
    }

    const label newI = fieldNames.size();

    fieldNames.setSize(newI + 1);
    fieldNames[newI] = field.name();

    auto& fieldPtrs = fields<Type>();
    fieldPtrs.setSize(newI + 1);
    fieldPtrs[newI] = &field;

    auto& fieldResiduals = residuals<Type>();
    fieldResiduals.setSize(newI + 1);
    allocateResidual(fieldResiduals[newI]);
}


template<class Type>
inline void dgMassProjection::addResidual
(
    const word& fieldName,
    const label cellID,
    const List<Type>& residual
)
{
    validateCellID(cellID);

    List<Type>& cellResidual = residualBuffer<Type>(fieldName)[cellID];

    if (cellResidual.size() != residual.size())
    {
        FatalErrorInFunction
            << "Cannot add residual to field " << fieldName
            << " in cell " << cellID << ": target has "
            << cellResidual.size() << " dofs but contribution has "
            << residual.size() << '.'
            << abort(FatalError);
    }

    forAll(cellResidual, dofI)
    {
        cellResidual[dofI] += residual[dofI];
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
