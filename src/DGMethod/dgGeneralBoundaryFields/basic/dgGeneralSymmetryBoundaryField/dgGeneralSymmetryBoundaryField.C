/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DGFoam: Discontinuous Galerkin CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | GPU-friendly CFD solver framework
     \\/     M anipulation  |
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
    along with DGFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dgGeneralSymmetryBoundaryField.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "tensor.H"
#include "VectorSpace.H"
#include "dgGeneralBoundaryFieldMacros.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeDgGeneralBoundaryField(dgGeneralSymmetryBoundaryField);

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * //

template<class Type>
dgGeneralSymmetryBoundaryField<Type>::dgGeneralSymmetryBoundaryField
(
    const word& name,
    const dictionary& dict
)
:
    dgGeneralBoundaryField<Type>(name, dict)
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<>
void dgGeneralSymmetryBoundaryField<scalar>::updateValue
(
    const label gaussID,
    const vector& n,
    const scalar& minusValue,
    const scalar& minusGrad,
    scalar& plusValue,
    scalar& plusGrad
) const
{
    plusValue = minusValue;
}

template<>
void dgGeneralSymmetryBoundaryField<vector>::updateValue
(
    const label gaussID,
    const vector& n,
    const vector& minusValue,
    const vector& minusGrad,
    vector& plusValue,
    vector& plusGrad
) const
{
    plusValue = minusValue - 2 * (n & minusValue) * n;
}

template<>
void dgGeneralSymmetryBoundaryField<tensor>::updateValue
(
    const label gaussID,
    const vector& n,
    const tensor& minusValue,
    const tensor& minusGrad,
    tensor& plusValue,
    tensor& plusGrad
) const
{
    tensor R = tensor::I - 2.0 * (n * n); // Reflection tensor, n * n is outer product of n with itself
    plusValue = R & minusValue & R;
}

template<>
void dgGeneralSymmetryBoundaryField<scalar>::updateGrad
(
    const label gaussID,
    const vector& n,
    const scalar& minusValue,
    const scalar& minusGrad,
    scalar& plusValue,
    scalar& plusGrad
) const
{
    plusGrad = minusGrad;
}

template<>
void dgGeneralSymmetryBoundaryField<vector>::updateGrad
(
    const label gaussID,
    const vector& n,
    const vector& minusValue,
    const vector& minusGrad,
    vector& plusValue,
    vector& plusGrad
) const
{
    plusGrad = minusGrad - 2 * (n & minusGrad) * n;
}

template<>
void dgGeneralSymmetryBoundaryField<tensor>::updateGrad
(
    const label gaussID,
    const vector& n,
    const tensor& minusValue,
    const tensor& minusGrad,
    tensor& plusValue,
    tensor& plusGrad
) const
{
    tensor R = tensor::I - 2.0 * (n * n); // Reflection tensor, n * n is outer product of n with itself
    plusGrad = R & minusGrad & R;
}

} // End namespace Foam

