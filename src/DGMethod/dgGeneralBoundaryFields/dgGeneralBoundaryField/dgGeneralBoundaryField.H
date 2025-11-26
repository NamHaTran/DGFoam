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
    Foam::dgGeneralBoundaryField

Description
    Base class for runtime-selectable DG boundary conditions.
    When creating a new boundary condition, please remember to override
    checkPatchType() method. Below are the PolyPatch types in OpenFOAM

    Physical boundary patches:
    - polyPatch               : Standard physical patch (generic)
    - wallPolyPatch           : Solid wall patch
    - emptyPolyPatch          : 2D/1D degenerate patch (no normal direction)
    - wedgePolyPatch          : Axisymmetric wedge patch
    - symmetryPolyPatch       : Legacy symmetry patch
    - symmetryPlanePolyPatch  : Modern symmetry patch
    
    Cyclic / periodic patches:
    - cyclicPolyPatch         : Simple periodic patch
    - cyclicAMIPolyPatch      : AMI periodic/rotating interface
    
    Region-coupling / mapping patches:
    - mappedPolyPatch         : Generic mapped patch (1:1 mapping)
    - mappedWallPolyPatch     : Mapped wall interface (CHT)
    - regionCoupledPolyPatch  : Patch coupling two regions (matching faces)
    - wallRegionCoupledPolyPatch : Region-coupled version of wall
    
    Parallel communication patch:
    - processorPolyPatch      : MPI processor patch
    
    Mesh-motion and utility patches (less commonly used):
    - faceZonePolyPatch       : Patch used for faceZone extrusions
    - ggiPolyPatch (older FOAM) : GGI interface (predecessor to cyclicAMI)

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgGeneralBoundaryField_H
#define Foam_dgGeneralBoundaryField_H

#include "word.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dgGeomMesh.H"
#include "fvPatch.H"

namespace Foam
{

// Forward declaration
template<class Type> class faceGaussField;

/*---------------------------------------------------------------------------*\
                    Class dgGeneralBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dgGeneralBoundaryField
{
protected:

    //- Reference to fvPatch (geometry of patch)
    const fvPatch& patch_;

    //- Reference to DG geometric mesh (for DG-level cell/face connectivity)
    const dgGeomMesh& dgMesh_;

    //- Patch dictionary used to read BC parameters
    const dictionary dict_;

public:

    //- Runtime type information
    TypeName("dgGeneralBoundaryField");

    //- Runtime selection table (dictionary constructor)
    declareRunTimeSelectionTable
    (
        autoPtr,
        dgGeneralBoundaryField,
        dictionary,
        (const fvPatch& patch, const dgGeomMesh& dgMesh, const dictionary& dict),
        (patch, dgMesh, dict)
    );

    // Constructors

        //- Construct from fvPatch, dgGeomMesh, and dictionary
        dgGeneralBoundaryField
        (
            const fvPatch& patch,
            const dgGeomMesh& dgMesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~dgGeneralBoundaryField() = default;

    // Factory Method

        static autoPtr<dgGeneralBoundaryField<Type>> New
        (
            const fvPatch& patch,
            const dgGeomMesh& dgMesh,
            const dictionary& dict
        );

    // Core Interface

        //- Update '+' side value at a DG cell boundary for Type scalar
        virtual void updateGhostStateValue
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const scalar& minusValue,
            scalar& plusValue
        ) const
        {
            plusValue = minusValue;
        }

        //- Update '+' side value at a DG cell boundary for Type vector
        virtual void updateGhostStateValue
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const vector& minusValue,
            vector& plusValue
        ) const
        {
            plusValue = minusValue;
        }

        //- Update '+' side value at a DG cell boundary for Type tensor
        virtual void updateGhostStateValue
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const tensor& minusValue,
            tensor& plusValue
        ) const
        {
            plusValue = minusValue;
        }

        //- Update '+' side gradient for Type scalar
        virtual void updateGhostStateGrad
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const scalar& minusGrad,
            scalar& plusGrad
        ) const
        {
            plusGrad = minusGrad;
        }

        //- Update '+' side gradient for Type vector
        virtual void updateGhostStateGrad
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const vector& minusGrad,
            vector& plusGrad
        ) const
        {
            plusGrad = minusGrad;
        }

        //- Update '+' side gradient for Type tensor
        virtual void updateGhostStateGrad
        (
            const label cellID,
            const label faceLocalID,
            const label localGauss,
            const vector& n,
            const tensor& minusGrad,
            tensor& plusGrad
        ) const
        {
            plusGrad = minusGrad;
        }

    // Access Functions

        inline const fvPatch& patch() const { return patch_; }

        inline const dgGeomMesh& dgMesh() const { return dgMesh_; }

        inline const dictionary& dict() const { return dict_; }

        inline const word& name() const { return patch_.name(); }

        inline word BCtype() const
        {
            word t;
            dict_.lookup("type") >> t;
            return t;
        }

    // Check type
        
        //- Verify patch type compatibility
        virtual void checkPatchType() const {}
};

} // End namespace Foam

/*---------------------------------------------------------------------------*\
                  Include template implementation
\*---------------------------------------------------------------------------*/

#include "dgGeneralBoundaryFieldI.H"

#endif
// ************************************************************************* //
