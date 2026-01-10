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
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::faceGaussField

Description
    Gauss-point field values on faces of a DG cell.

    Stores values at Gauss integration points on all faces of a cell,
    including both minus (owner) and plus (neighbour) sides.

SourceFiles
    faceGaussField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceGaussField_H
#define Foam_faceGaussField_H

#include "dgGeomMesh.H"
#include "dgGeomCell.H"
#include "dgGeomFace.H"
#include "cellDof.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceGaussField Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::faceGaussField
 * \brief Gauss-point field on faces of a DG cell
 *
 * \tparam Type Field value type (e.g. scalar, vector, tensor)
 *
 * \details
 * faceGaussField stores field values evaluated at Gauss integration
 * points on all faces of a DG cell.
 *
 * For each Gauss point on a face, the class stores:
 * - Minus-side value (owner cell)
 * - Plus-side value (neighbour cell or boundary)
 *
 * The class is designed for cell-centric DG formulations and
 * supports flux computation across cell interfaces.
 *
 * This class does not own geometry or DoF storage; it holds
 * non-owning pointers to dgGeomMesh, dgGeomCell, dgGeomFace,
 * and cellDof<Type>.
 */
template<class Type>
class faceGaussField
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** Pointer to DG geometric mesh (non-owning) */
        const dgGeomMesh* mesh_;

        /** Cell index associated with this face field */
        label cellID_;

        /** Pointer to DG geometric cell (non-owning) */
        const dgGeomCell* cell_;

        /** Pointer to list of face IDs belonging to the cell (non-owning) */
        const labelList* facesID_;

        /** Number of faces of the cell */
        label nFaces_;

        /** Number of Gauss points per face */
        label nGaussPerFace_;

        /** Total number of Gauss points on all faces */
        label nGauss_;

        /** List of DG geometric faces (non-owning) */
        List<dgGeomFace*> faces_;

        /** Pointers to DoF containers of owner and neighbour cells */
        List<const cellDof<Type>*> cellsDof_;

        /** Offset of Gauss points for each face */
        List<label> gaussOffset_;

        /** Owner patch IDs for each face */
        List<label> ownerPatchIDs_;

        /** Face normal vectors at Gauss points */
        List<vector> n_;

        /** Face areas */
        List<scalar> faceAreas_;

        /** Values on plus (neighbour) side */
        List<Type> plusValues_;

        /** Values on minus (owner) side */
        List<Type> minusValues_;

        /** Flux vectors at Gauss points */
        List<vector> flux_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /** Default constructor */
        faceGaussField();

        /**
         * \brief Construct from DoF pointers and mesh
         *
         * Initializes face Gauss-point storage using DoFs
         * from owner and neighbour cells.
         */
        faceGaussField
        (
            List<const cellDof<Type>*>& cellsDof,
            const dgGeomMesh* mesh
        );

        /** Copy constructor */
        faceGaussField(const faceGaussField<Type>& other);

        /**
         * \brief Construct with uniform initial value
         */
        faceGaussField
        (
            const label cellID,
            const dgGeomMesh* mesh,
            const Type& initVal
        );

        /**
         * \brief Construct with zero-initialized values
         */
        faceGaussField
        (
            const label cellID,
            const dgGeomMesh* mesh
        );

    //====================================================================//
    //                              Access
    //====================================================================//

        /** Return cell index */
        inline label cellID() const { return cellID_; }

        /**
         * \brief Return DG geometric mesh
         *
         * Fatal error if pointer is null.
         */
        inline const dgGeomMesh* dgMesh() const
        {
            if (!mesh_)
            {
                FatalErrorInFunction << "mesh_ pointer is null" << abort(FatalError);
            }
            return mesh_;
        }

        /** Return number of faces */
        inline label nFaces() const { return nFaces_; }

        /** Return total number of Gauss points */
        inline label nGauss() const { return nGauss_; }

        /** Return number of Gauss points per face */
        inline label nGaussPerFace() const { return nGaussPerFace_; }

        /** Return list of DG geometric faces */
        inline const List<dgGeomFace*>& faces() const
        {
            return faces_;
        }

    //====================================================================//
    //                      Gauss-point value access
    //====================================================================//

        /** Access plus-side value by global Gauss index */
        inline Type& plusValueAt(const label gi)       { return plusValues_[gi]; }

        /** Access minus-side value by global Gauss index */
        inline Type& minusValueAt(const label gi)      { return minusValues_[gi]; }

        /** Const access to plus-side value */
        inline const Type& plusValue(const label gi) const  { return plusValues_[gi]; }

        /** Const access to minus-side value */
        inline const Type& minusValue(const label gi) const { return minusValues_[gi]; }

        /** Access plus-side value by local face and Gauss index */
        inline Type& plusValueOnFace(const label fI, const label gI)
        {
            return plusValues_[gaussOffset_[fI] + gI];
        }

        inline const Type& plusValueOnFace(const label fI, const label gI) const
        {
            return plusValues_[gaussOffset_[fI] + gI];
        }

        /** Access minus-side value by local face and Gauss index */
        inline Type& minusValueOnFace(const label fI, const label gI)
        {
            return minusValues_[gaussOffset_[fI] + gI];
        }

        inline const Type& minusValueOnFace(const label fI, const label gI) const
        {
            return minusValues_[gaussOffset_[fI] + gI];
        }

        /** Access minus-side value using global face ID */
        inline Type& minusValueOnGlobalFace(const label faceID, const label gI)
        {
            return minusValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        inline const Type& minusValueOnGlobalFace(const label faceID, const label gI) const
        {
            return minusValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        /** Access plus-side value using global face ID */
        inline Type& plusValueOnGlobalFace(const label faceID, const label gI)
        {
            return plusValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        inline const Type& plusValueOnGlobalFace(const label faceID, const label gI) const
        {
            return plusValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        /** Return minus-side values */
        inline const List<Type>& minusSide() const { return minusValues_; }
        inline List<Type>& minusSide() { return minusValues_; }

        /** Return plus-side values */
        inline const List<Type>& plusSide() const { return plusValues_; }
        inline List<Type>& plusSide() { return plusValues_; }

        /** Access flux on a face Gauss point */
        inline const vector& fluxOnFace(const label fI, const label gI) const
        {
            return flux_[gaussOffset_[fI] + gI];
        }

        inline vector& fluxOnFace(const label fI, const label gI)
        {
            return flux_[gaussOffset_[fI] + gI];
        }

    //====================================================================//
    //                      Geometry access
    //====================================================================//

        /** Return face normals */
        inline const List<vector>& normals() const { return n_; }

        /** Return face areas */
        inline const List<scalar>& faceAreas() const { return faceAreas_; }

        /** Return global face ID from local face index */
        inline label globalFaceID(const label fI) const { return (*facesID_)[fI]; }

        /**
         * \brief Return local face index from global face ID
         *
         * Fatal error if face is not found.
         */
        inline label localFaceIndex(const label faceID) const 
        {
            for (label fI = 0; fI < nFaces_; ++fI)
            {
                if ((*facesID_)[fI] == faceID)
                {
                    return fI;
                }
            }

            FatalErrorInFunction
                << "Face ID " << faceID << " not found in cell " << cellID_
                << abort(FatalError);
            return -1;
        }

        /** Check whether face is boundary */
        inline bool isBoundary(const label fI) const { return faces_[fI]->isBoundary(); }

        /** Check whether input cell is owner of face */
        inline bool isOwner(const label fI, const label cellID) const 
        { 
            return faces_[fI]->isOwner(cellID); 
        }

        /** Return owner patch ID */
        inline label getOwnerPatchID(const label fI) const
        {
            return ownerPatchIDs_[fI];
        }

        /** Set DoF pointers */
        inline void setCellsDof(const List<const cellDof<Type>*>& lst)
        {
            cellsDof_ = lst;
        }

    //====================================================================//
    //                        Interpolation
    //====================================================================//

        /**
         * \brief Interpolate DoF values to face Gauss points
         */
        void interpolateFromDof();

    //====================================================================//
    //                           Operators
    //====================================================================//

        /** Copy assignment */
        faceGaussField<Type>& operator=(const faceGaussField<Type>& other);

        /** Assign uniform value */
        faceGaussField<Type>& operator=(const Type& v);

        /** Initialize flux storage */
        inline void initFlux() 
        { 
            if (flux_.size() == 0)
                flux_.setSize(nGauss_); 
        }

    //====================================================================//
    //                             Output
    //====================================================================//

        /** Stream output of face Gauss values */
        friend Ostream& operator<<
        (
            Ostream& os,
            const faceGaussField<Type>& fgf
        )
        {
            os << "faceGaussField for cell " << fgf.cellID() << nl;

            for (label i = 0; i < fgf.nGauss(); ++i)
            {
                os << "  Gauss pt " << i << " minus: " << fgf.minusValue(i) << nl;
                os << "  Gauss pt " << i << " plus : " << fgf.plusValue(i)  << nl;
            }
            return os;
        }

    //====================================================================//
    //                Component access (vector fields)
    //====================================================================//

        /** Return x-component */
        faceGaussField<scalar> x() const;

        /** Return y-component */
        faceGaussField<scalar> y() const;

        /** Return z-component */
        faceGaussField<scalar> z() const;

        /** Set x-component */
        void setX(const faceGaussField<scalar>& fx);

        /** Set y-component */
        void setY(const faceGaussField<scalar>& fy);

        /** Set z-component */
        void setZ(const faceGaussField<scalar>& fz);

        /** Set all components */
        void setXYZ
        (
            const faceGaussField<scalar>& fx,
            const faceGaussField<scalar>& fy,
            const faceGaussField<scalar>& fz
        );

    //====================================================================//
    //                              Limiting
    //====================================================================//

        /** Limit maximum value */
        void limitMax(const scalar maxVal);

        /** Limit minimum value */
        void limitMin(const scalar minVal);

        /** Limit value range */
        void limitRange(const scalar minVal, const scalar maxVal);
};

} // End namespace Foam

#include "faceGaussFieldI.H"
#include "faceGaussFieldScalarI.H"
#include "faceGaussFieldVectorI.H"

#endif

// ************************************************************************* //
