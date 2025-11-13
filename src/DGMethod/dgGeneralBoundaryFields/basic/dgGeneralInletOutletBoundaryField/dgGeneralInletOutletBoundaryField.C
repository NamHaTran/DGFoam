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

#include "dgGeneralInletOutletBoundaryField.H"
#include "dgGeneralBoundaryFieldMacros.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Instantiation * * * * * * * * * * * * * * //

makeDgGeneralBoundaryField(dgGeneralInletOutletBoundaryField);

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

template<class Type>
dgGeneralInletOutletBoundaryField<Type>::dgGeneralInletOutletBoundaryField
(
    const word& name,
    const dictionary& dict
)
:
    dgGeneralBoundaryField<Type>(name, dict)
{
    dict.readEntry("inletValue", inletValue_);
}

// * * * * * * * * * * * * * * updateValue * * * * * * * * * * * * * * * * //

template<class Type>
void dgGeneralInletOutletBoundaryField<Type>::updateValue
(
    const label gaussID,
    const vector& n,
    const Type& minusValue,
    const Type& minusGrad,
    Type& plusValue,
    Type& plusGrad
) const
{
    /*
    const faceGaussField<vector>& V =
        this->ctxPtr_->template lookupFaceField<vector>("U");

    const vector& Vn = V.plusValue(gaussID);
    scalar dotProd = Vn & n;

    if (dotProd < 0)  // Inlet: flow into domain
    {
        plusValue = 2 * inletValue_ - minusValue;
    }
    else // Outlet: flow out of domain
    {
        plusValue = minusValue;
    */
}


// * * * * * * * * * * * * * * updateGrad * * * * * * * * * * * * * * * * * //

template<class Type>
void dgGeneralInletOutletBoundaryField<Type>::updateGrad
(
    const label gaussID,
    const vector& n,
    const Type& minusValue,
    const Type& minusGrad,
    Type& plusValue,
    Type& plusGrad
) const
{
    plusGrad = minusGrad;
}


// * * * * * * * * * * * * * * End * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
