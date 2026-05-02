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

Class
    Foam::faceGaussField

Description
    Gauss-point field values on faces of a DG cell.

    Stores values at Gauss integration points on all faces of a cell,
    including both minus (cell-local interior) and plus (exterior) sides.

SourceFiles
    faceGaussFieldI.H
    faceGaussFieldScalarI.H
    faceGaussFieldVectorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceGaussField_H
#define Foam_faceGaussField_H

#include "dgGeomMesh.H"
#include "dgGeomCell.H"
#include "dgGeomFace.H"
#include "cellDof.H"
#include "List.H"
#include "boundaryGaussField.H"
#include "refCount.H"

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
 * - Minus-side value (the current cell / interior trace)
 * - Plus-side value (the neighbour or boundary / exterior trace)
 * - Middle-state value used by weak face formulations
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
:
    public refCount
{
private:

    //====================================================================//
    //                           Data members
    //====================================================================//

        /** \brief Non-owning pointer to the DG geometric mesh. */
        const dgGeomMesh* mesh_;

        /** \brief Mesh cell index associated with this face Gauss field. */
        label cellID_;

        /** \brief Non-owning pointer to the DG geometric cell. */
        const dgGeomCell* cell_;

        /** \brief Non-owning pointer to the global face IDs of the cell. */
        const labelList* facesID_;

        /** \brief Number of faces belonging to the cell. */
        label nFaces_;

        /** \brief Number of Gauss points stored per face. */
        label nGaussPerFace_;

        /** \brief Total number of stored face Gauss points. */
        label nGauss_;

        /** \brief Non-owning pointers to the DG geometric faces of the cell. */
        List<dgGeomFace*> faces_;

        /** \brief Pointers to the owner and neighbour cell DoF containers. */
        List<const cellDof<Type>*> cellsDof_;

        /** \brief Local copies of neighbour DoFs received on processor faces. */
        List<cellDof<Type>> bufferedNbrDof_;

        /** \brief Validity flags for the buffered processor-neighbour DoFs. */
        List<bool> bufferedNbrDofValid_;

        /** \brief Flat-list offset of the first Gauss point of each face. */
        List<label> gaussOffset_;

        /** \brief Owner patch IDs for the faces, or -1 for internal faces. */
        List<label> ownerPatchIDs_;

        /** \brief Outward face normals with respect to the owner cell. */
        List<vector> n_;

        /** \brief Face areas. */
        List<scalar> faceAreas_;

        /** \brief Exterior-trace values at face Gauss points. */
        List<Type> plusValues_;

        /** \brief Interior-trace values at face Gauss points. */
        List<Type> minusValues_;

        /** \brief Middle-state values at face Gauss points. */
        List<Type> midValues_;

        /** \brief Flux vectors stored at face Gauss points when requested. */
        List<vector> flux_;

public:

    //====================================================================//
    //                           Constructors
    //====================================================================//

        /** \brief Construct an empty face Gauss field. */
        faceGaussField();

        /**
         * \brief Construct from DoF pointers and mesh
         *
         * Initializes face Gauss-point storage using DoFs
         * from owner and neighbour cells.
         *
         * \param cellsDof Owner-cell and neighbour-cell DoF pointers.
         * \param mesh DG geometric mesh containing the cell and its faces.
         */
        faceGaussField
        (
            List<const cellDof<Type>*>& cellsDof,
            const dgGeomMesh* mesh
        );

        /**
         * \brief Copy constructor.
         *
         * \param other Source face Gauss field to copy from.
         */
        faceGaussField(const faceGaussField<Type>& other);

        /**
         * \brief Move constructor.
         *
         * \param other Source face Gauss field whose storage is moved from.
         */
        faceGaussField(faceGaussField<Type>&& other) noexcept;

        /**
         * \brief Construct with a uniform initial value.
         *
         * \param cellID Mesh cell index.
         * \param mesh DG geometric mesh containing the cell.
         * \param initVal Initial value assigned at every face Gauss point.
         */
        faceGaussField
        (
            const label cellID,
            const dgGeomMesh* mesh,
            const Type& initVal
        );

        /**
         * \brief Construct with zero-initialized values.
         *
         * \param cellID Mesh cell index.
         * \param mesh DG geometric mesh containing the cell.
         */
        faceGaussField
        (
            const label cellID,
            const dgGeomMesh* mesh
        );

        /**
         * \brief Construct by consuming a temporary faceGaussField.
         *
         * \param t Temporary face Gauss field whose storage is consumed.
         */
        faceGaussField(const tmp<faceGaussField<Type>>& t);

    //====================================================================//
    //                              Access
    //====================================================================//

        /** \brief Return the associated mesh cell index. */
        inline label cellID() const { return cellID_; }

        /**
         * \brief Return DG geometric mesh
         *
         * Fatal error if pointer is null.
         *
         * \return Pointer to the DG geometric mesh.
         */
        inline const dgGeomMesh* dgMesh() const
        {
            if (!mesh_)
            {
                FatalErrorInFunction << "mesh_ pointer is null" << abort(FatalError);
            }
            return mesh_;
        }

        /** \brief Return the number of faces in the cell. */
        inline label nFaces() const { return nFaces_; }

        /** \brief Return the total number of stored face Gauss points. */
        inline label nGauss() const { return nGauss_; }

        /** \brief Return the number of stored Gauss points per face. */
        inline label nGaussPerFace() const { return nGaussPerFace_; }

        /** \brief Return the list of DG geometric faces of the cell. */
        inline const List<dgGeomFace*>& faces() const
        {
            return faces_;
        }

    //====================================================================//
    //                      Gauss-point value access
    //====================================================================//

        /** \brief Return writable access to one flat plus-side value. */
        inline Type& plusValueAt(const label gi)       { return plusValues_[gi]; }

        /** \brief Return writable access to one flat minus-side value. */
        inline Type& minusValueAt(const label gi)      { return minusValues_[gi]; }

        /** \brief Return read-only access to one flat plus-side value. */
        inline const Type& plusValue(const label gi) const  { return plusValues_[gi]; }

        /** \brief Return read-only access to one flat minus-side value. */
        inline const Type& minusValue(const label gi) const { return minusValues_[gi]; }

        /** \brief Return writable access to one flat middle-state value. */
        inline Type& midValueAt(const label gi)       { return midValues_[gi]; }

        /** \brief Return read-only access to one flat middle-state value. */
        inline const Type& midValue(const label gi) const  { return midValues_[gi]; }

        /** \brief Return writable access to a plus-side value on one face. */
        inline Type& plusValueOnFace(const label fI, const label gI)
        {
            return plusValues_[gaussOffset_[fI] + gI];
        }

        /** \brief Return read-only access to a plus-side value on one face. */
        inline const Type& plusValueOnFace(const label fI, const label gI) const
        {
            return plusValues_[gaussOffset_[fI] + gI];
        }

        /** \brief Return writable access to a minus-side value on one face. */
        inline Type& minusValueOnFace(const label fI, const label gI)
        {
            return minusValues_[gaussOffset_[fI] + gI];
        }

        /** \brief Return read-only access to a minus-side value on one face. */
        inline const Type& minusValueOnFace(const label fI, const label gI) const
        {
            return minusValues_[gaussOffset_[fI] + gI];
        }

        /** \brief Return writable access to a middle-state value on one face. */
        inline Type& midValueOnFace(const label fI, const label gI)
        {
            return midValues_[gaussOffset_[fI] + gI];
        }

        /** \brief Return read-only access to a middle-state value on one face. */
        inline const Type& midValueOnFace(const label fI, const label gI) const
        {
            return midValues_[gaussOffset_[fI] + gI];
        }

        /** \brief Return writable access to a minus-side value using a global face ID. */
        inline Type& minusValueOnGlobalFace(const label faceID, const label gI)
        {
            return minusValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        /** \brief Return read-only access to a minus-side value using a global face ID. */
        inline const Type& minusValueOnGlobalFace(const label faceID, const label gI) const
        {
            return minusValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        /** \brief Return writable access to a plus-side value using a global face ID. */
        inline Type& plusValueOnGlobalFace(const label faceID, const label gI)
        {
            return plusValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        /** \brief Return read-only access to a plus-side value using a global face ID. */
        inline const Type& plusValueOnGlobalFace(const label faceID, const label gI) const
        {
            return plusValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        /** \brief Return writable access to a middle-state value using a global face ID. */
        inline Type& midValueOnGlobalFace(const label faceID, const label gI)
        {
            return midValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        /** \brief Return read-only access to a middle-state value using a global face ID. */
        inline const Type& midValueOnGlobalFace(const label faceID, const label gI) const
        {
            return midValues_[gaussOffset_[localFaceIndex(faceID)] + gI];
        }

        /** \brief Return read-only access to the full minus-side storage. */
        inline const List<Type>& minusSide() const { return minusValues_; }
        /** \brief Return writable access to the full minus-side storage. */
        inline List<Type>& minusSide() { return minusValues_; }

        /** \brief Return read-only access to the full plus-side storage. */
        inline const List<Type>& plusSide() const { return plusValues_; }
        /** \brief Return writable access to the full plus-side storage. */
        inline List<Type>& plusSide() { return plusValues_; }

        /** \brief Return read-only access to the full middle-state storage. */
        inline const List<Type>& midValues() const { return midValues_; }
        /** \brief Return writable access to the full middle-state storage. */
        inline List<Type>& midValues() { return midValues_; }

        /** \brief Return read-only access to one flux vector on a face. */
        inline const vector& fluxOnFace(const label fI, const label gI) const
        {
            return flux_[gaussOffset_[fI] + gI];
        }

        /** \brief Return writable access to one flux vector on a face. */
        inline vector& fluxOnFace(const label fI, const label gI)
        {
            return flux_[gaussOffset_[fI] + gI];
        }

    //====================================================================//
    //                      Geometry access
    //====================================================================//

        /** \brief Return the outward face normals of the owner cell. */
        inline const List<vector>& normals() const { return n_; }

        /** \brief Return the face areas. */
        inline const List<scalar>& faceAreas() const { return faceAreas_; }

        /** \brief Return the global face ID from a local face index. */
        inline label globalFaceID(const label fI) const { return (*facesID_)[fI]; }

        /**
         * \brief Return local face index from global face ID
         *
         * Fatal error if face is not found.
         *
         * \param faceID Global face ID.
         * \return Local face index inside the owner cell.
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

        /** \brief Return true if the specified local face is a boundary face. */
        inline bool isBoundary(const label fI) const { return faces_[fI]->isBoundary(); }

        /** \brief Return true if the specified local face lies on a processor patch. */
        inline bool isProcessorPatch(const label fI) const
        {
            return faces_[fI]->isProcessorPatch();
        }

        /** \brief Return true if the input cell is the owner of the specified face. */
        inline bool isOwner(const label fI, const label cellID) const 
        { 
            return faces_[fI]->isOwner(cellID); 
        }

        /** \brief Return the owner patch ID, or -1 for internal faces. */
        inline label getOwnerPatchID(const label fI) const
        {
            return ownerPatchIDs_[fI];
        }

        /** \brief Replace the owner/neighbour DoF pointer list. */
        inline void setCellsDof(const List<const cellDof<Type>*>& lst)
        {
            cellsDof_ = lst;
        }

        /**
         * \brief Interpolate only owner-side local faces from DoFs.
         *
         * Local internal faces for which this cell is the neighbour are left
         * untouched so a field-level pass can copy the already interpolated
         * owner trace instead of reconstructing the same face twice.
         */
        void interpolateOwnerFacesFromDof();

        /**
         * \brief Fill this neighbour-side local face from its owner cell face.
         *
         * The owner plus trace is this cell's minus trace, and the owner minus
         * trace is this cell's plus trace at the same physical Gauss points.
         */
        void copyNeighbourFaceFromOwner
        (
            const label localFaceI,
            const faceGaussField<Type>& ownerFaceField
        );

        /**
         * \brief Store neighbour DoFs received for one processor face.
         *
         * The input is copied into face-local owned storage so callers can
         * safely access the data after the communication helper returns.
         *
         * \param localFaceI Local face index.
         * \param nbrDof Neighbour-cell DoFs received through synchronization.
         */
        inline void setBufferedNeighbourDof
        (
            const label localFaceI,
            const cellDof<Type>& nbrDof
        )
        {
            bufferedNbrDof_[localFaceI] = nbrDof;
            bufferedNbrDofValid_[localFaceI] = true;
        }

        /**
         * \brief Access neighbour DoFs for one face.
         *
         * Internal faces return the directly linked neighbour cell DoFs.
         * Processor faces return the last communicated neighbour DoFs.
         *
         * \param localFaceI Local face index.
         * \return Read-only neighbour cell DoF data.
         */
        inline const cellDof<Type>& neighbourCellDof
        (
            const label localFaceI
        ) const
        {
            if (localFaceI < 0 || localFaceI >= nFaces_)
            {
                FatalErrorInFunction
                    << "Face index " << localFaceI << " is out of range [0, "
                    << nFaces_ - 1 << "] for cell " << cellID_
                    << abort(FatalError);
            }

            const dgGeomFace& face = *faces_[localFaceI];

            if (!face.isBoundary())
            {
                const cellDof<Type>* nbrDofPtr = cellsDof_[localFaceI + 1];

                if (!nbrDofPtr)
                {
                    FatalErrorInFunction
                        << "Missing neighbour DoF pointer on internal face "
                        << localFaceI << " of cell " << cellID_
                        << abort(FatalError);
                }

                return *nbrDofPtr;
            }

            if (!face.isProcessorPatch())
            {
                FatalErrorInFunction
                    << "Face " << localFaceI << " of cell " << cellID_
                    << " is a physical boundary face and has no neighbour DoF."
                    << abort(FatalError);
            }

            if (!bufferedNbrDofValid_[localFaceI])
            {
                FatalErrorInFunction
                    << "Processor neighbour DoF buffer is not initialized for "
                    << "face " << localFaceI << " of cell " << cellID_
                    << ". Synchronize processor boundaries before accessing "
                    << "neighbourCellDof()."
                    << abort(FatalError);
            }

            return bufferedNbrDof_[localFaceI];
        }

        /**
         * \brief Extract ghost-state values on boundary faces into a flat
         * boundaryGaussField using local-face order then Gauss-point order.
         */
        boundaryGaussField<Type> extractBCGhostState() const;

        /**
         * \brief Assign flat boundaryGaussField values back to ghost-state
         * values on boundary faces using the same ordering as
         * extractBCGhostState().
         */
        void assignBCGhostState(const boundaryGaussField<Type>& bgf);

        /**
         * \brief Extract middle-state values on boundary faces into a flat
         * boundaryGaussField using local-face order then Gauss-point order.
         */
        boundaryGaussField<Type> extractMidValue() const;

        /**
         * \brief Assign flat boundaryGaussField values back to middle-state
         * values on boundary faces using the same ordering as
         * extractMidValue().
         */
        void assignMidValue(const boundaryGaussField<Type>& bgf);

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

        /** \brief Copy assignment operator. */
        faceGaussField<Type>& operator=(const faceGaussField<Type>& other);

        /** \brief Move assignment operator. */
        faceGaussField<Type>& operator=
        (
            faceGaussField<Type>&& other
        ) noexcept;

        /** \brief Assign by consuming a temporary faceGaussField. */
        faceGaussField<Type>& operator=
        (
            const tmp<faceGaussField<Type>>& t
        );

        /** \brief Assign a uniform value to all plus, minus, and middle Gauss points. */
        faceGaussField<Type>& operator=(const Type& v);

        /** \brief Allocate the flux storage if it has not been allocated yet. */
        inline void initFlux() 
        { 
            if (flux_.size() == 0)
                flux_.setSize(nGauss_); 
        }

    //====================================================================//
    //                             Output
    //====================================================================//

        /**
         * \brief Write the face Gauss-point values to an output stream.
         *
         * \param os Output stream.
         * \param fgf Face Gauss field to serialize.
         * \return Reference to the output stream.
         */
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

        /** \brief Return the x-component as a scalar faceGaussField. */
        faceGaussField<scalar> x() const;

        /** \brief Return the y-component as a scalar faceGaussField. */
        faceGaussField<scalar> y() const;

        /** \brief Return the z-component as a scalar faceGaussField. */
        faceGaussField<scalar> z() const;

        /** \brief Replace the x-component using a scalar faceGaussField. */
        void setX(const faceGaussField<scalar>& fx);

        /** \brief Replace the y-component using a scalar faceGaussField. */
        void setY(const faceGaussField<scalar>& fy);

        /** \brief Replace the z-component using a scalar faceGaussField. */
        void setZ(const faceGaussField<scalar>& fz);

        /** \brief Replace all vector components using scalar faceGaussFields. */
        void setXYZ
        (
            const faceGaussField<scalar>& fx,
            const faceGaussField<scalar>& fy,
            const faceGaussField<scalar>& fz
        );

    //====================================================================//
    //                              Limiting
    //====================================================================//

        /** \brief Limit plus/minus values to a maximum value. */
        void limitMax(const scalar maxVal);

        /** \brief Limit plus/minus values to a minimum value. */
        void limitMin(const scalar minVal);

        /** \brief Limit plus/minus values to the range [minVal, maxVal]. */
        void limitRange(const scalar minVal, const scalar maxVal);
};

} // End namespace Foam

#include "faceGaussFieldI.H"
#include "faceGaussFieldScalarI.H"
#include "faceGaussFieldVectorI.H"

#endif

// ************************************************************************* //
