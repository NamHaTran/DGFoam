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
      serves as an intermediate coordinate system between eta-space
      and the Montjoie reference space.
    - Montjoie reference coordinate system (x-space):
      used for the evaluation of basis functions.

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
     * \brief Transform reference coordinates from xi-space to x-space for hexahedral elements
     * \details Montjoie reference space for hexahedral element is within [0, 1] in all directions
     * Reference: https://www.math.u-bordeaux.fr/~durufle/montjoie/hexahedra.php
     * 
     *
     * \param xiPt Gauss point in reference space (xi)
     */
    inline Foam::vector xiToXHex(const Foam::vector& xiPt)
    {
        Foam::vector xPt = xiPt;
        xPt.x() = 0.5 * (xiPt.x() + 1.0);
        xPt.y() = 0.5 * (xiPt.y() + 1.0);
        xPt.z() = 0.5 * (xiPt.z() + 1.0);

        return xPt;
    };

    /**
     * \brief Transform reference coordinates from x-space to xi-space for hexahedral elements
     * \details Inverse mapping of the Montjoie reference hexahedral element.
     *          Montjoie reference space is within [0, 1] in all directions.
     *
     * \param xPt Point in Montjoie reference space (x)
     */
    inline Foam::vector xToXiHex(const Foam::vector& xPt)
    {
        Foam::vector xiPt = xPt;
        xiPt.x() = 2.0 * xPt.x() - 1.0;
        xiPt.y() = 2.0 * xPt.y() - 1.0;
        xiPt.z() = 2.0 * xPt.z() - 1.0;

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
     * \brief Transform reference coordinates from xi-space to x-space for prismatic elements
     * \details Montjoie reference space for prismatic element is within [0, 1] in all directions
     * Reference: https://www.math.u-bordeaux.fr/~durufle/montjoie/wedge.php
     * 
     *
     * \param xiPt Gauss point in reference space (xi)
     */
    inline Foam::vector xiToXPrism(const Foam::vector& xiPt)
    {
        Foam::scalar xi1 = xiPt.x();
        Foam::scalar xi2 = xiPt.y();
        Foam::scalar xi3 = xiPt.z();

        Foam::vector xPt = xiPt;
        xPt.x() = (1 + xi1)*(1 - xi3)/4;
        xPt.y() = 0.5 * (xi3 + 1.0);
        xPt.z() = 0.5 * (xi2 + 1.0);

        return xPt;
    };

    /**
     * \brief Transform reference coordinates from x-space to xi-space for prismatic elements
     * \details Inverse mapping of the Montjoie reference prismatic element.
     *          Montjoie reference space is within [0, 1] in all directions.
     *
     * \param xPt Point in Montjoie reference space (x)
     */
    inline Foam::vector xToXiPrism(const Foam::vector& xPt)
    {
        const Foam::scalar x = xPt.x();
        const Foam::scalar y = xPt.y();
        const Foam::scalar z = xPt.z();

        Foam::vector xiPt = xPt;

        // xi3 from y
        xiPt.z() = 2.0 * y - 1.0;

        // xi2 from z
        xiPt.y() = 2.0 * z - 1.0;

        // xi1 from x and y (using 1 - xi3 = 2(1 - y))
        xiPt.x() = 2.0 * x / (1.0 - y) - 1.0;

        return xiPt;
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
     * \brief Transform reference coordinates from xi-space to x-space for pyramidic elements
     * \details Montjoie reference space for pyramidic element is within [0, 1] in all directions
     * Reference: https://www.math.u-bordeaux.fr/~durufle/montjoie/pyramid.php
     * 
     *
     * \param xiPt Gauss point in reference space (xi)
     */
    inline Foam::vector xiToXPyramid(const Foam::vector& xiPt)
    {
        Foam::scalar xi1 = xiPt.x();
        Foam::scalar xi2 = xiPt.y();
        Foam::scalar xi3 = xiPt.z();

        Foam::vector xPt = xiPt;
        xPt.x() = (1 - xi3) * xi1/2;
        xPt.y() = (1 - xi3) * xi2/2;
        xPt.z() = 0.5 * (xi3 + 1.0);

        return xPt;
    };

    /**
     * \brief Transform reference coordinates from x-space to xi-space for pyramidic elements
     * \details Inverse mapping of the Montjoie reference pyramidic element.
     *          Montjoie reference space is within [0, 1] in all directions.
     *
     * \param xPt Point in Montjoie reference space (x)
     */
    inline Foam::vector xToXiPyramid(const Foam::vector& xPt)
    {
        const Foam::scalar x = xPt.x();
        const Foam::scalar y = xPt.y();
        const Foam::scalar z = xPt.z();

        Foam::vector xiPt = xPt;

        xiPt.z() = 2.0 * y - 1.0;
        xiPt.y() = y / (1 - z);
        xiPt.x() = x / (1 - z);

        return xiPt;
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

    /**
     * \brief Transform reference coordinates from xi-space to x-space for tetrahedral elements
     * \details Montjoie reference space for tetrahedral element is within [0, 1] in all directions
     * Reference: https://www.math.u-bordeaux.fr/~durufle/montjoie/tetrahedra.php
     * 
     *
     * \param xiPt Gauss point in reference space (xi)
     */
    inline Foam::vector xiToXTet(const Foam::vector& xiPt)
    {
        return xiToXHex(xiPt);
    };

    /**
     * \brief Transform reference coordinates from x-space to xi-space for tetrahedral elements
     * \details Inverse mapping of the Montjoie reference tetrahedral element.
     *          Montjoie reference space is within [0, 1] in all directions.
     *
     * \param xPt Point in Montjoie reference space (x)
     */
    inline Foam::vector xToXiTet(const Foam::vector& xPt)
    {
        return xToXiHex(xPt);
    };

} // End namespace Foam

#endif

// ************************************************************************* //
