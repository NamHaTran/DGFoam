/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2024-2026 Ha Nam Tran
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
    Foam::TVBMinmod

Description
    TVB-minmod troubled-cell detector based on face-average jumps.

\*---------------------------------------------------------------------------*/

#ifndef Foam_TVBMinmod_H
#define Foam_TVBMinmod_H

#include "dgField.H"
#include "troubleCellDetector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class TVBMinmod Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::TVBMinmod
 * \brief Face-based TVB-minmod detector for conservative DG fields.
 *
 * \details
 * This detector compares the owner-side average on each face with the cell
 * average and limits that jump by a two-argument TVB-minmod operator.
 * For one scalar-like quantity, the candidate jump is
 *
 * \f[
 *     a_0 = \bar{u}_{f,-} - \bar{u}_K,
 * \f]
 *
 * where \f$\bar{u}_{f,-}\f$ is the owner-side face average and
 * \f$\bar{u}_K\f$ is the cell average.
 *
 * For an internal face, the comparison jump is taken from the neighboring
 * cell average,
 *
 * \f[
 *     a_1 = \bar{u}_{K_{nb}} - \bar{u}_K.
 * \f]
 *
 * For a boundary face, the exterior state reconstructed on that face may be
 * used instead when \c includeBoundaryFaces is enabled.
 *
 * The relaxed TVB-minmod check is
 *
 * \f[
 *     \tilde{m}(a_0, a_1; Mh^2),
 * \f]
 *
 * with cell size \f$h\f$ and user-specified TVB constant \f$M\f$. A cell is
 * marked troubled when the relaxed value differs from the owner-face jump.
 *
 * The detector loops over the fields listed in \c checkFields.
 * Scalar fields are checked directly, while vector fields are checked
 * component-by-component. The cell is marked as troubled as soon as one
 * field component violates the TVB-minmod check.
 *
 * The following optional dictionary entries are supported:
 * \verbatim
 * checkFields          (rho U E);
 * M                    1.0;
 * includeBoundaryFaces true;
 * report               false;
 * \endverbatim
 */
class TVBMinmod
:
    public troubleCellDetector
{
    // Private data

        /** TVB constant appearing in the relaxation threshold \f$Mh^2\f$. */
    scalar M_;

        /** Include boundary faces by comparing owner and exterior face states. */
    bool includeBoundaryFaces_;

public:

    /** \brief Runtime type information. */
    TypeName("TVBMinmod");

    /**
     * \brief Construct the TVB-minmod detector.
     *
     * \param dict Detector dictionary
     * \param mesh DG mesh providing geometry and field registry access
     */
    TVBMinmod
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /**
     * \brief Evaluate the TVB-minmod sensor for one cell.
     *
     * \param cellID Target cell index
     * \return True if any conservative component violates the TVB-minmod check
     */
    bool detect(const label cellID) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
