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

#include "dgField.H"
#include "dgCellType.H"
#include "troubleCellDetector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PerssonPeraire Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::PerssonPeraire
 * \brief Modal-decay detector following Persson and Peraire.
 *
 * \details
 * This detector's idea:
 * the highest-order modal contribution of a checked field component, meaning the complete
 * modal subspace added when the approximation order increases from
 * \f$p-1\f$ to \f$p\f$, is compared against the total modal energy of the
 * reconstructed field component inside one cell.
 *
 * In the present DG framework the detector obtains:
 * - modal coefficients from the checked field \c dofField,
 * - basis values, quadrature weights, and Jacobians from \c dgGeomCell.
 *
 * The resulting smoothness sensor is
 *
 * \f[
 *     S_e = \frac{\int_K (u_p)^2\,dV}{\int_K u^2\,dV},
 * \f]
 *
 * where \f$u_p\f$ is the reconstructed contribution of all modes that
 * belong to the order-\f$p\f$ shell for the selected scalar quantity or
 * vector component.
 *
 * The detector then evaluates
 *
 * \f[
 *     s_e = \log_{10}(S_e),
 * \f]
 *
 * and marks a cell once \f$s_e\f$ enters the activation band
 * \f$s_e \ge s_0 - \kappa\f$.
 *
 * The following optional dictionary entries are supported:
 * \verbatim
 * checkFields   (rho);
 * s0           <scalar>;   // default: log10(1/(p+1)^4)
 * kappa        1.0;
 * report       false;
 * \endverbatim
 */
class PerssonPeraire
:
    public troubleCellDetector
{
    // Private data

        /** Expected number of modal coefficients for tetrahedra. */
        label tetNDof_;

        /** Expected number of modal coefficients for hexahedra. */
        label hexNDof_;

        /** Expected number of modal coefficients for prisms. */
        label prismNDof_;

        /** Expected number of modal coefficients for pyramids. */
        label pyramidNDof_;

        /** Cached order-\f$p\f$ modal shell indices for tetrahedra. */
        List<label> tetPModeDof_;

        /** Cached order-\f$p\f$ modal shell indices for hexahedra. */
        List<label> hexPModeDof_;

        /** Cached order-\f$p\f$ modal shell indices for prisms. */
        List<label> prismPModeDof_;

        /** Cached order-\f$p\f$ modal shell indices for pyramids. */
        List<label> pyramidPModeDof_;

        /** Logarithmic activation center of the Persson-Peraire sensor. */
        scalar s0_;

        /** Half-width of the smooth activation band in log-sensor space. */
        scalar kappa_;

        /** \brief Return the cached order-\f$p\f$ modal shell indices for a cell type. */
        const List<label>& pModeDof(const dgCellType cellType) const;

        /** \brief Return the expected modal coefficient count for a cell type. */
        label expectedNDof(const dgCellType cellType) const;

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
 * \return True if any checked scalar field or vector component enters the configured band
 */
    virtual bool detect(const label cellID) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
