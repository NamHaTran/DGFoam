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
    Foam::KXRCF

Description
    KXRCF troubled-cell detector based on inflow-face jumps.

\*---------------------------------------------------------------------------*/

#ifndef Foam_KXRCF_H
#define Foam_KXRCF_H

#include "dgField.H"
#include "troubleCellDetector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class KXRCF Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::KXRCF
 * \brief Inflow-face jump detector following the KXRCF smoothness sensor.
 *
 * \details
 * This detector approximates the KXRCF indicator
 *
 * \f[
 *     I_j = \int_{\partial \Omega_j^-}(Q_j - Q_{nb})\,ds,
 * \f]
 *
 * where \f$\partial \Omega_j^-\f$ is the inflow part of the cell boundary,
 * identified from the sign of the local normal direction indicator
 * \f$\boldsymbol{d}\cdot \boldsymbol{n}\f$.
 *
 * The normalized sensor is evaluated as
 *
 * \f[
 *     T_j =
 *     \frac{|I_j|}
 *          {h^{(p+1)/2}\,|\partial \Omega_j^-|\,\|Q_j\|_\infty},
 * \f]
 *
 * and the cell is marked troubled when \f$T_j > T_0\f$.
 *
 * Optionally, the detector can switch to the LPR-style normalization
 * proposed by Li et al.,
 *
 * \f[
 *     S_j =
 *     \frac{|I_j|}
 *          {h\,|\partial \Omega_j^-|\,|\bar{Q}_j|},
 * \f]
 *
 * where \f$\bar{Q}_j\f$ is the cell-average value of the monitored quantity.
 *
 * The threshold must be a positive scalar. Smaller values make the detector
 * more aggressive because more cells satisfy the troubled-cell criterion,
 * whereas larger values make it more selective.
 *
 * Practical guidance:
 * - For the original KXRCF normalization (\c LPRMode false), \c threshold=1
 *   is the canonical choice from the standard sensor definition.
 * - For the LPR normalization (\c LPRMode true), Li et al. use a constant of
 *   order one. A value near \c 1 is a reasonable default for simplex-like
 *   elements, while larger values such as \c 2 may be more appropriate when
 *   following quadrilateral/tensor-product calibrations.
 *
 * In short: choose \c threshold > 0, start from \c 1.0, decrease it to catch
 * weaker oscillations, and increase it to avoid over-limiting smooth regions.
 *
 * The following optional dictionary entries are supported:
 * \verbatim
 * checkFields          (rho);
 * directionField       rhoU;
 * LPRMode              false;
 * threshold            1.0;
 * includeBoundaryFaces true;
 * report               false;
 * \endverbatim
 *
 * Header note on \c directionField:
 * - \c directionField is any DG vector field used only to decide whether a
 *   face quadrature point is inflow or outflow through the sign of
 *   \f$\boldsymbol{d}\cdot\boldsymbol{n}\f$.
 * - The detector does not normalize the sensor with this field and does not
 *   measure jumps of this field unless it is also listed separately in
 *   \c checkFields.
 * - In a compressible conservative solver, \c rhoU is often a good choice
 *   because \f$(\rho\boldsymbol{U})\cdot\boldsymbol{n}\f$ has the same sign as
 *   \f$\boldsymbol{U}\cdot\boldsymbol{n}\f$ when \f$\rho > 0\f$, while avoiding
 *   an extra primitive-field refresh before limiting.
 * - If a case already maintains a consistent primitive velocity field before
 *   detection, \c U can also be used as the \c directionField.
 */
class KXRCF
:
    public troubleCellDetector
{
    // Private data

        /**
         * Vector field used only to determine the sign of the face-normal
         * direction indicator when classifying inflow/outflow.
         */
        const dgField<vector>& directionField_;

        /** Switch between original KXRCF and the Li-et-al. LPR normalization. */
        bool LPRMode_;

        /** Positive cutoff: smaller is more aggressive, larger is more selective. */
        scalar threshold_;

        /** Include physical boundary faces in the inflow-face integral. */
        bool includeBoundaryFaces_;

public:

    /** \brief Runtime type information. */
    TypeName("KXRCF");

    /**
     * \brief Construct the KXRCF detector.
     *
     * \param dict Detector dictionary
     * \param mesh DG mesh providing geometry and field registry access
     */
    KXRCF
    (
        const dictionary& dict,
        const dgGeomMesh& mesh
    );

    /**
 * \brief Evaluate the KXRCF sensor for one cell.
 *
 * \param cellID Target cell index
 * \return True if any checked scalar field or vector component exceeds threshold
 */
    virtual bool detect(const label cellID) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
