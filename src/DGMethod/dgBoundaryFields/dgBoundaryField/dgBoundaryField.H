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
    Foam::dgBoundaryField

Description
    Base class for runtime-selectable DG boundary conditions.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgBoundaryField_H
#define Foam_dgBoundaryField_H

#include "word.H"
#include "dictionary.H"
#include "vector.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fieldsContext.H"

namespace Foam
{

//- Boundary-condition enforcement form
enum class dgBCForm : unsigned char
{
    strong = 0,   // Strong imposition
    weak   = 1    // Weak/penalty imposition
};


/*---------------------------------------------------------------------------*\
                           Class dgBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgBoundaryField
{
protected:

    // Attributes

        //- Name of boundary patch
        word name_;

        //- Dictionary for this boundary
        dictionary dict_;

        //- Pointer to the fields context for registration/lookup
        const fieldsContext* ctxPtr_ = nullptr;

        //- Enforcement form (strong/weak)
        dgBCForm bcForm_;


    // Helpers

        //- Parse "form" from dictionary (default = strong)
        static dgBCForm parseForm_(const dictionary& dict);


public:

    //- Runtime type info
    TypeName("dgBoundaryField");

    //- Runtime selection support
    declareRunTimeSelectionTable
    (
        autoPtr,
        dgBoundaryField,
        dictionary,
        (const word& name, const dictionary& dict),
        (name, dict)
    );

    // Constructors

        //- Construct from name and dictionary
        dgBoundaryField(const word& name, const dictionary& dict);

    //- Destructor
    virtual ~dgBoundaryField() = default;

    // Factory

        //- Static factory method
        static autoPtr<dgBoundaryField<Type>> New
        (
            const word& name,
            const dictionary& dict
        );

    // Core interface

        //- Update '+' side value according to BC (may use form)
        virtual void updateValue
        (
            const vector& n,
            const Type& minusValue,
            const Type& minusGrad,
            Type& plusValue,
            Type& plusGrad
        ) const = 0;

        //- Update '+' side gradient according to BC (may use form)
        virtual void updateGrad
        (
            const vector& n,
            const Type& minusValue,
            const Type& minusGrad,
            Type& plusValue,
            Type& plusGrad
        ) const = 0;

    // Context

        //- Set the fields context
        inline void setContext(const fieldsContext& ctx) { ctxPtr_ = &ctx; }

        //- Get the fields context pointer (may be null)
        inline const fieldsContext* context() const { return ctxPtr_; }

    // Form accessors

        //- Return enforcement form
        inline dgBCForm form() const { return bcForm_; }

        //- True if strong form
        inline bool isStrong() const { return bcForm_ == dgBCForm::strong; }

        //- True if weak form
        inline bool isWeak()   const { return bcForm_ == dgBCForm::weak; }

        //- Return form as word ("strong" or "weak")
        word formName() const;
};

} // End namespace Foam

//#ifdef NoRepository
//    #include "dgBoundaryField.C"
//#endif

#endif

// ************************************************************************* //
