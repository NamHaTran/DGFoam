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
    Foam::modifiedMinMod

Description
    Modified minmod troubled-cell detector adapted from the legacy 2D code.

\*---------------------------------------------------------------------------*/

#ifndef Foam_modifiedMinMod_H
#define Foam_modifiedMinMod_H

#include "troubleCellDetector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class modifiedMinMod Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::modifiedMinMod
 * \brief Face-based troubled-cell detector derived from the legacy p-adaptive limiter.
 *
 * \details
 * The legacy 2D code used a modified-minmod check as part of the p-adaptive
 * limiting logic. In that algorithm, one compares a directional jump measured
 * from the current cell polynomial against the neighboring-cell means and
 * applies a TVB-like relaxation
 *
 * \f[
 *     \tilde{m}(a_1,\dots,a_n; M h^2 + \epsilon).
 * \f]
 *
 * In the current 3D DG code the same idea is adapted to the data structures
 * available in \c dgGeomCell and \c GaussField:
 * - cell means are taken directly from the zeroth conservative modal
 *   coefficient, which represents the cell average in the modal DG basis,
 * - face states are approximated by averaging owner-side face Gauss values,
 * - the detector is applied independently to density, each momentum component,
 *   and total energy.
 *
 * If any conservative component violates the modified-minmod check, the cell
 * is marked as troubled.
 */
class modifiedMinMod
:
    public troubleCellDetector
{
public:

    /** \brief Runtime type information. */
    TypeName("modifiedMinMod");

    /**
     * \brief Construct the modified-minmod detector.
     *
     * \param dict Detector dictionary
     * \param mesh DG mesh providing geometry and field registry access
     */
    modifiedMinMod
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /**
     * \brief Evaluate the modified-minmod troubled-cell check for one cell.
     *
     * \param cellID Target cell index
     * \return True if any conservative component triggers the detector
     */
    virtual bool detect(const label cellID) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
