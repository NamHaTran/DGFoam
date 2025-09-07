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

Class
    Foam::dgInletOutletBoundaryField

Description
    General inletOutlet boundary condition for DG fields.

    - Requires a velocity field named "U" (of type faceGaussField<vector>)
      to be registered in the fieldsContext.
    - Checks flow direction using U·n at each Gauss point.
    - If flow is incoming (U·n < 0), applies a fixedValue condition:
        \f[
            \phi^+ = 2 \cdot \phi_{\text{inlet}} - \phi^-
        \f]
    - If flow is outgoing (U·n >= 0), applies a zeroGradient condition:
        \f[
            \phi^+ = \phi^-
        \f]
    - Gradient is always set as:
        \f[
            \nabla \phi^+ = \nabla \phi^-
        \f]

    - The value "inletValue" must be specified in the dictionary.

SourceFiles
    dgInletOutletBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef dgGeneralInletOutletBoundaryField_H
#define dgGeneralInletOutletBoundaryField_H

#include "dgGeneralBoundaryField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dgGeneralInletOutletBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgGeneralInletOutletBoundaryField
:
    public dgGeneralBoundaryField<Type>
{
protected:

    //- Inlet value for Dirichlet BC
    Type inletValue_;

public:

    //- Runtime type information
    TypeName("generalInletOutlet");

    //- Constructor from name and dictionary
    dgGeneralInletOutletBoundaryField
    (
        const word& name,
        const dictionary& dict
    );
    
    //- Update value at plus side
    virtual void updateValue
    (
        const label gaussID,
        const vector& n,
        const Type& minusValue,
        const Type& minusGrad,
        Type& plusValue,
        Type& plusGrad
    ) const override;

    //- Update gradient at plus side
    virtual void updateGrad
    (
        const label gaussID,
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