/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2025 Ha Nam Tran
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

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::dgSourceTerm<Type>::dgSourceTerm
(
    const word& termName,
    const dgBasisField& basis,
    const cellGaussField<Type>& S
)
:
    termName_(termName),
    cellID_(basis.cellID()),
    basis_(basis),
    S_(new cellGaussField<Type>(S)),
    nDof_(basis_.nDof()),
    R_(nDof_, pTraits<Type>::zero)
{}


template<class Type>
Foam::dgSourceTerm<Type>::dgSourceTerm
(
    const word& termName,
    const dgBasisField& basis,
    const tmp<cellGaussField<Type>>& tS
)
:
    termName_(termName),
    cellID_(basis.cellID()),
    basis_(basis),
    S_(tS),
    nDof_(basis_.nDof()),
    R_(nDof_, pTraits<Type>::zero)
{}


template<class Type>
Foam::dgSourceTerm<Type>::dgSourceTerm
(
    const word& termName,
    const dgBasisField& basis,
    const GaussField<Type>& S
)
:
    termName_(termName),
    cellID_(basis.cellID()),
    basis_(basis),
    S_(new cellGaussField<Type>(S.cellField())),
    nDof_(basis_.nDof()),
    R_(nDof_, pTraits<Type>::zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::dgSourceTerm<Type>::calcResidual()
{
    forAll (R_, i)
    {
        R_[i] = pTraits<Type>::zero;
    }

    const cellGaussField<Type>& S = S_();

    if (S.cellID() != cellID_)
    {
        FatalErrorInFunction
            << "Source field cellID " << S.cellID()
            << " does not match basis cellID " << cellID_
            << " for source term " << termName_
            << abort(FatalError);
    }

    const List<scalar>& w = S.cell()->weights();
    const List<scalar>& J3D = S.cell()->J3D();

    for (label dofI = 0; dofI < nDof_; ++dofI)
    {
        const cellGaussField<scalar>& phi =
            basis_.getBasis(dofI).cellField();

        Type volTerm = pTraits<Type>::zero;

        for (label gpI = 0; gpI < S.size(); ++gpI)
        {
            volTerm += (J3D[gpI]*w[gpI]*phi[gpI])*S[gpI];
        }

        R_[dofI] += volTerm;
    }
}

// ************************************************************************* //
