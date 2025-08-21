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

  Class dgGeneralSymmetryBoundaryField Declaration

  Description
    - Implements a general symmetry boundary condition for scalar, vector, and
      tensor fields.

    - For scalar fields:
        * Value at the + side is equal to the - side.
        * Gradient is preserved.

    - For vector fields:
        * Value is reflected across the plane normal to the boundary normal:
              U+ = U- - 2(n·U-)n
        * Gradient is also reflected accordingly:
              S+ = S- - 2(n·S-)n

    - For tensor fields:
        * Value and gradient are both reflected using:
              T+ = R·T-·R
              R = I - 2(n⊗n)

    - This condition assumes the boundary is a symmetry plane:
        * No flux passes through the boundary.
        * Useful for reducing domain size in symmetric problems.

  Usage
    boundaryField
    {
        patchName
        {
            type    generalSymmetry;
        }
    }

  Note
    - Only geometric information (normal vector n) and internal values are
      used. No input parameters required from dictionary.

SourceFiles
    dgGeneralZeroGradientBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralSymmetryBoundaryField_H
#define dgGeneralSymmetryBoundaryField_H

#include "dgGeneralBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dgGeneralSymmetryBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgGeneralSymmetryBoundaryField
:
    public dgGeneralBoundaryField<Type>
{
public:

    //- Runtime type information
    TypeName("generalSymmetry");

    //- Constructor from name and dictionary
    dgGeneralSymmetryBoundaryField
    (
        const word& name,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgGeneralSymmetryBoundaryField() = default;

    //- Update '+' value according to symmetry condition
    virtual void updateValue
    (
        const vector& n,
        const Type& minusValue,
        const Type& minusGrad,
        Type& plusValue,
        Type& plusGrad
    ) const override;

    //- Update '+' gradient according to symmetry condition
    virtual void updateGrad
    (
        const vector& n,
        const Type& minusValue,
        const Type& minusGrad,
        Type& plusValue,
        Type& plusGrad
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //

