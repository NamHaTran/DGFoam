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
    Foam::PerssonPeraire

Description
    Persson-Peraire troubled-cell detector based on density modal decay.

\*---------------------------------------------------------------------------*/

#ifndef Foam_PerssonPeraire_H
#define Foam_PerssonPeraire_H

#include "troubleCellDetector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PerssonPeraire Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::PerssonPeraire
 * \brief Density-based modal-decay detector following Persson and Peraire.
 *
 * \details
 * This detector's idea:
 * the highest-order modal contribution of density is compared against the
 * total modal energy of the reconstructed density field inside one cell.
 *
 * In the present DG framework the detector obtains:
 * - modal coefficients from the density \c dofField,
 * - basis values, quadrature weights, and Jacobians from \c dgGeomCell.
 *
 * The resulting smoothness sensor is
 *
 * \f[
 *     S_e = \frac{\int_K (\rho_p)^2\,dV}{\int_K \rho^2\,dV},
 * \f]
 *
 * where \f$\rho_p\f$ is the contribution of the highest modal coefficient.
 * A cell is flagged when \f$S_e\f$ exceeds a polynomial-order-dependent
 * threshold.
 */
class PerssonPeraire
:
    public troubleCellDetector
{
public:

    /** \brief Runtime type information. */
    TypeName("PerssonPeraire");

    /**
     * \brief Construct the Persson-Peraire detector.
     *
     * \param dict Detector dictionary
     * \param mesh DG mesh providing geometry and field registry access
     */
    PerssonPeraire
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /**
     * \brief Evaluate the Persson-Peraire sensor for one cell.
     *
     * \param cellID Target cell index
     * \return True if the smoothness sensor exceeds the configured threshold
     */
    virtual bool detect(const label cellID) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
