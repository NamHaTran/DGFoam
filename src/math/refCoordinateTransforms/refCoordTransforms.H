/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

Description
    A collection of utilities for transforming between different reference
    coordinate systems.

    DGFoam uses three reference coordinate systems:
    - DG hexahedral reference coordinate system (eta-space):
      used for numerical integration.
    - DG element reference coordinate system (xi-space):
      used for shape function evaluation and interpolation.

\*---------------------------------------------------------------------------*/

#ifndef Foam_refCoordTransforms_H
#define Foam_refCoordTransforms_H

#include "scalar.H"
#include "vector.H"

namespace Foam
{
    /**
     * \brief Transform reference coordinates from eta-space to xi-space for hexahedral elements
     * \details No transformation
     * 
     *
     * \param etaPt Gauss point in reference space (eta)
     */
    inline Foam::vector etaToXiHex(const Foam::vector& etaPt)
    {
        return etaPt;
    };

    /**
     * \brief Transform reference coordinates from xi-space to eta-space for hexahedral elements
     * \details No transformation
     *
     * \param xiPt Point in reference space (xi)
     */
    inline Foam::vector xiToEtaHex(const Foam::vector& xiPt)
    {
        return xiPt;
    };

    /**
     * \brief Transform reference coordinates from eta-space to xi-space for prismatic elements
     * \details collapse the eta1 direction
     * 
     *
     * \param etaPt Gauss point in reference space (eta)
     */
    inline Foam::vector etaToXiPrism(const Foam::vector& etaPt)
    {
        Foam::scalar eta1 = etaPt.x();
        Foam::scalar eta3 = etaPt.z();

        Foam::vector xiPt = etaPt;
        xiPt.x() = (1 + eta1)*(1 - eta3)/2 - 1.0;

        return xiPt;
    };

    /**
     * \brief Transform reference coordinates from xi-space to eta-space for prismatic elements
     * \details collapse the eta1 direction
     * 
     *
     * \param xiPt Gauss point in reference space (xi)
     */
    inline Foam::vector xiToEtaPrism(const Foam::vector& xiPt)
    {
        Foam::scalar xi1 = xiPt.x();
        Foam::scalar xi3 = xiPt.z();

        Foam::vector etaPt = xiPt;
        etaPt.x() = 2.0 * (xi1 + 1.0)/(1 - xi3) - 1.0;

        return etaPt;
    };

    /**
     * \brief Transform reference coordinates from eta-space to xi-space for pyramidic elements
     * \details collapse the eta1 and eta2 direction
     * 
     *
     * \param etaPt Gauss point in reference space (eta)
     */
    inline Foam::vector etaToXiPyramid(const Foam::vector& etaPt)
    {
        Foam::scalar eta1 = etaPt.x();
        Foam::scalar eta2 = etaPt.y();
        Foam::scalar eta3 = etaPt.z();

        Foam::vector xiPt = etaPt;
        xiPt.x() = (1 + eta1)*(1 - eta3)/2 - 1.0;
        xiPt.y() = (1 + eta2)*(1 - eta3)/2 - 1.0;

        return xiPt;
    };

    /**
     * \brief Transform reference coordinates from xi-space to eta-space for pyramidic elements
     * \details collapse the eta1 and eta2 direction
     * 
     *
     * \param xiPt Gauss point in reference space (xi)
     */
    inline Foam::vector xiToEtaPyramid(const Foam::vector& xiPt)
    {
        Foam::scalar xi1 = xiPt.x();
        Foam::scalar xi2 = xiPt.y();
        Foam::scalar xi3 = xiPt.z();

        Foam::vector etaPt = xiPt;
        etaPt.x() = 2.0 * (xi1 + 1.0)/(1 - xi3) - 1.0;
        etaPt.y() = 2.0 * (xi2 + 1.0)/(1 - xi3) - 1.0;

        return etaPt;
    };

    /**
     * \brief Transform reference coordinates from eta-space to xi-space for tetrahedral elements
     * \details collapse the eta1 and eta2 direction
     * 
     *
     * \param etaPt Gauss point in reference space (eta)
     */
    inline Foam::vector etaToXiTet(const Foam::vector& etaPt)
    {
        const Foam::scalar eta1 = etaPt.x();
        const Foam::scalar eta2 = etaPt.y();
        const Foam::scalar eta3 = etaPt.z();

        Foam::vector xiPt = etaPt;

        // xi3 = eta3
        xiPt.z() = eta3;

        // xi2 from eta2 and eta3
        // eta2 + 1 = 2*(xi2 + 1)/(1 - xi3)
        xiPt.y() = 0.5 * (eta2 + 1.0) * (1.0 - eta3) - 1.0;

        // xi1 from eta1, eta2, eta3
        // eta1 + 1 = 2*(xi1 + 1)/(-xi2 - xi3)
        // with: -xi2 - xi3 = (1 - eta3)(1 - eta2)/2
        xiPt.x() =
            0.25 * (eta1 + 1.0) * (1.0 - eta3) * (1.0 - eta2) - 1.0;

        return xiPt;
    };

    /**
     * \brief Transform reference coordinates from xi-space to eta-space for tetrahedral elements
     * \details collapse the eta1 and eta2 direction
     * 
     *
     * \param xiPt Gauss point in reference space (xi)
     */
    inline Foam::vector xiToEtaTet(const Foam::vector& xiPt)
    {
        Foam::scalar xi1 = xiPt.x();
        Foam::scalar xi2 = xiPt.y();
        Foam::scalar xi3 = xiPt.z();

        Foam::vector etaPt = xiPt;
        etaPt.x() = 2.0 * (xi1 + 1.0)/(-xi2 - xi3) - 1.0;
        etaPt.y() = 2.0 * (xi2 + 1.0)/(1 - xi3) - 1.0;

        return etaPt;
    };

} // End namespace Foam

#endif

// ************************************************************************* //
