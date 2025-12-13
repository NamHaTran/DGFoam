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
#include "symmetryPolyPatch.H"
#include "vector.H"
#include "tensor.H"
#include "VectorSpace.H"

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
        const fvPatch& patch,
        const dgGeomMesh& dgMesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dgGeneralSymmetryBoundaryField() = default;

    //- Member methods
        inline void updateGhostStateValue
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const scalar& minusValue,
            scalar& plusValue
        ) const override
        {
            plusValue = minusValue;
        }

        void updateGhostStateValue
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const vector& minusValue,
            vector& plusValue
        ) const override
        {
            plusValue = minusValue - 2 * (n & minusValue) * n;
        }

        /*
        void dgGeneralSymmetryBoundaryField<tensor>::updateGhostStateValue
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const vector& minusValue,
            vector& plusValue
        ) const override
        {
            tensor R = tensor::I - 2.0 * n * n; // Reflection tensor
            plusValue = R & minusValue & R;
        }

        void dgGeneralSymmetryBoundaryField<vector>::updateGhostStateGrad
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const vector& minusGrad,
            vector& plusGrad
        ) const override
        {
            plusGrad = minusGrad - 2 * (n & minusGrad) * n;
        }

        void dgGeneralSymmetryBoundaryField<tensor>::updateGhostStateGrad
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const vector& minusGrad,
            vector& plusGrad
        ) const override
        {
            tensor R = tensor::I - 2.0 * n * n; // Reflection tensor
            plusGrad = R & minusGrad & R;
        }
        */

        //- Check patch type
        void checkPatchType() const override
        {
            if (this->patch_.type() != "symmetry")
            {
                FatalErrorInFunction
                    << "Boundary condition " << type() <<" can only "
                    << "be applied to patch type: \n"
                    << "    " << "symmetry \n"
                    << "but patch "
                    << this->patch_.name() << " is of type " << this->patch_.type() << nl
                    << exit(FatalError);
            }
        }
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgGeneralSymmetryBoundaryFieldI.H"

#endif

// ************************************************************************* //

