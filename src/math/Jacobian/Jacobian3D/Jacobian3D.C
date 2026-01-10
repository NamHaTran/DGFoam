/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "Jacobian3D.H"
#include "error.H"
#include <cmath>
#include "scalar.H"
#include "label.H"
#include "List.H"

namespace Foam
{
namespace math
{

Foam::scalar Foam::HexJacobian3D
(
    const vector& gaussPt,             // Gauss point in reference space (eta)
    const List<vector>& cellVertices   // Physical coordinates of hex vertices
)
{
    // ---------------------------------------------------------------------
    // Extract reference coordinates (eta-space)
    // ---------------------------------------------------------------------
    const scalar eta1 = gaussPt.x();
    const scalar eta2 = gaussPt.y();
    const scalar eta3 = gaussPt.z();

    // ---------------------------------------------------------------------
    // Gradients of shape functions w.r.t eta (dPhi/dEta)
    // Ordering: A B C D E F G H
    // ---------------------------------------------------------------------
    List<vector> gradPhi_Eta(cellVertices.size(), Zero);

    // Vertex A
    gradPhi_Eta[0].x() = -(1 - eta2)*(1 - eta3)/8;
    gradPhi_Eta[0].y() = -(1 - eta1)*(1 - eta3)/8;
    gradPhi_Eta[0].z() = -(1 - eta1)*(1 - eta2)/8;

    // Vertex B
    gradPhi_Eta[1].x() =  (1 - eta2)*(1 - eta3)/8;
    gradPhi_Eta[1].y() = -(1 + eta1)*(1 - eta3)/8;
    gradPhi_Eta[1].z() = -(1 + eta1)*(1 - eta2)/8;

    // Vertex C
    gradPhi_Eta[2].x() = -(1 + eta2)*(1 - eta3)/8;
    gradPhi_Eta[2].y() =  (1 - eta1)*(1 - eta3)/8;
    gradPhi_Eta[2].z() = -(1 - eta1)*(1 + eta2)/8;

    // Vertex D
    gradPhi_Eta[3].x() =  (1 + eta2)*(1 - eta3)/8;
    gradPhi_Eta[3].y() =  (1 + eta1)*(1 - eta3)/8;
    gradPhi_Eta[3].z() = -(1 + eta1)*(1 + eta2)/8;

    // Vertex E
    gradPhi_Eta[4].x() = -(1 - eta2)*(1 + eta3)/8;
    gradPhi_Eta[4].y() = -(1 - eta1)*(1 + eta3)/8;
    gradPhi_Eta[4].z() =  (1 - eta1)*(1 - eta2)/8;

    // Vertex F
    gradPhi_Eta[5].x() =  (1 - eta2)*(1 + eta3)/8;
    gradPhi_Eta[5].y() = -(1 + eta1)*(1 + eta3)/8;
    gradPhi_Eta[5].z() =  (1 + eta1)*(1 - eta2)/8;

    // Vertex G
    gradPhi_Eta[6].x() = -(1 + eta2)*(1 + eta3)/8;
    gradPhi_Eta[6].y() =  (1 - eta1)*(1 + eta3)/8;
    gradPhi_Eta[6].z() =  (1 - eta1)*(1 + eta2)/8;

    // Vertex H
    gradPhi_Eta[7].x() =  (1 + eta2)*(1 + eta3)/8;
    gradPhi_Eta[7].y() =  (1 + eta1)*(1 + eta3)/8;
    gradPhi_Eta[7].z() =  (1 + eta1)*(1 + eta2)/8;

    // ---------------------------------------------------------------------
    // Compute physical gradients: d(x,y,z)/d(eta1,eta2,eta3)
    // Each vector stores partial derivatives w.r.t eta
    // ---------------------------------------------------------------------
    vector gradX1_Eta(Zero); // d(x)/d(eta)
    vector gradX2_Eta(Zero); // d(y)/d(eta)
    vector gradX3_Eta(Zero); // d(z)/d(eta)

    for (label i = 0; i < cellVertices.size(); ++i)
    {
        gradX1_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].x();
        gradX1_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].x();
        gradX1_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].x();

        gradX2_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].y();
        gradX2_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].y();
        gradX2_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].y();

        gradX3_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].z();
        gradX3_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].z();
        gradX3_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].z();
    }

    // ---------------------------------------------------------------------
    // Jacobian tensor: columns are gradients of x,y,z w.r.t eta
    // ---------------------------------------------------------------------
    const tensor J3D
    (
        gradX1_Eta,
        gradX2_Eta,
        gradX3_Eta
    );

    // ---------------------------------------------------------------------
    // Return determinant of Jacobian
    // ---------------------------------------------------------------------
    return det(J3D);
}

Foam::scalar Foam::PrismJacobian3D
(
    const vector& gaussPt,             // Gauss point in reference space (eta)
    const List<vector>& cellVertices   // Physical coordinates of prism vertices
)
{
    // ---------------------------------------------------------------------
    // Extract reference coordinates (eta-space)
    // ---------------------------------------------------------------------
    const scalar eta1 = gaussPt.x();
    const scalar eta2 = gaussPt.y();
    const scalar eta3 = gaussPt.z();

    // ---------------------------------------------------------------------
    // Gradients of shape functions w.r.t eta (dPhi/dEta)
    //
    // Vertex correspondence:
    //
    //   cellVertices id -> prism vertex -> gradPhi_Eta id
    //   -------------------------------------------------
    //        0        ->       E        ->       0
    //        1        ->       A        ->       1
    //        2        ->       B        ->       2
    //        3        ->       F        ->       3
    //        4        ->       D        ->       4
    //        5        ->       C        ->       5
    //
    // ---------------------------------------------------------------------
    List<vector> gradPhi_Eta(cellVertices.size(), Zero);

    // ---------------------------------------------------------------------
    // Vertex A (cellVertices[1])
    // ---------------------------------------------------------------------
    gradPhi_Eta[1].x() = -(1 - eta2)*(1 - eta3)/8;
    gradPhi_Eta[1].y() = -(1 - eta1)*(1 - eta3)/8;
    gradPhi_Eta[1].z() = -(1 - eta1)*(1 - eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex B (cellVertices[2])
    // ---------------------------------------------------------------------
    gradPhi_Eta[2].x() =  (1 - eta2)*(1 - eta3)/8;
    gradPhi_Eta[2].y() = -(1 + eta1)*(1 - eta3)/8;
    gradPhi_Eta[2].z() = -(1 + eta1)*(1 - eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex C (cellVertices[5])
    // ---------------------------------------------------------------------
    gradPhi_Eta[5].x() = -(1 + eta2)*(1 - eta3)/8;
    gradPhi_Eta[5].y() =  (1 - eta1)*(1 - eta3)/8;
    gradPhi_Eta[5].z() = -(1 - eta1)*(1 + eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex D (cellVertices[4])
    // ---------------------------------------------------------------------
    gradPhi_Eta[4].x() =  (1 + eta2)*(1 - eta3)/8;
    gradPhi_Eta[4].y() =  (1 + eta1)*(1 - eta3)/8;
    gradPhi_Eta[4].z() = -(1 + eta1)*(1 + eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex E (cellVertices[0])
    // ---------------------------------------------------------------------
    gradPhi_Eta[0].x() =  0;
    gradPhi_Eta[0].y() = -(1 + eta3)/4;
    gradPhi_Eta[0].z() =  (1 - eta2)/4;

    // ---------------------------------------------------------------------
    // Vertex F (cellVertices[3])
    // ---------------------------------------------------------------------
    gradPhi_Eta[3].x() =  0;
    gradPhi_Eta[3].y() =  (1 + eta3)/4;
    gradPhi_Eta[3].z() =  (1 + eta2)/4;

    // ---------------------------------------------------------------------
    // Compute physical gradients: d(x,y,z)/d(eta1,eta2,eta3)
    // ---------------------------------------------------------------------
    vector gradX1_Eta(Zero); // d(x)/d(eta)
    vector gradX2_Eta(Zero); // d(y)/d(eta)
    vector gradX3_Eta(Zero); // d(z)/d(eta)

    for (label i = 0; i < cellVertices.size(); ++i)
    {
        gradX1_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].x();
        gradX1_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].x();
        gradX1_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].x();

        gradX2_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].y();
        gradX2_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].y();
        gradX2_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].y();

        gradX3_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].z();
        gradX3_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].z();
        gradX3_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].z();
    }

    // ---------------------------------------------------------------------
    // Jacobian tensor (columns = gradients of x, y, z w.r.t eta)
    // ---------------------------------------------------------------------
    const tensor J3D
    (
        gradX1_Eta,
        gradX2_Eta,
        gradX3_Eta
    );

    // ---------------------------------------------------------------------
    // Return determinant of Jacobian
    // ---------------------------------------------------------------------
    return det(J3D);
}

Foam::scalar Foam::PyramidJacobian3D
(
    const vector& gaussPt,             // Gauss point in reference space (eta)
    const List<vector>& cellVertices   // Physical coordinates of pyramid vertices
)
{
    // ---------------------------------------------------------------------
    // Extract reference coordinates (eta-space)
    // ---------------------------------------------------------------------
    const scalar eta1 = gaussPt.x();
    const scalar eta2 = gaussPt.y();
    const scalar eta3 = gaussPt.z();

    // ---------------------------------------------------------------------
    // Gradients of shape functions w.r.t eta (dPhi/dEta)
    //
    // Vertex correspondence:
    //
    //   cellVertices id -> pyramid vertex -> gradPhi_Eta id
    //   ---------------------------------------------------
    //        0        ->        A         ->        0
    //        1        ->        B         ->        1
    //        2        ->        C         ->        2
    //        3        ->        D         ->        3
    //        4        ->        E         ->        4
    //
    // ---------------------------------------------------------------------
    List<vector> gradPhi_Eta(cellVertices.size(), Zero);

    // ---------------------------------------------------------------------
    // Vertex A (cellVertices[0])
    // ---------------------------------------------------------------------
    gradPhi_Eta[0].x() = -(1 - eta2)*(1 - eta3)/8;
    gradPhi_Eta[0].y() = -(1 - eta1)*(1 - eta3)/8;
    gradPhi_Eta[0].z() = -(1 - eta1)*(1 - eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex B (cellVertices[1])
    // ---------------------------------------------------------------------
    gradPhi_Eta[1].x() =  (1 - eta2)*(1 - eta3)/8;
    gradPhi_Eta[1].y() = -(1 + eta1)*(1 - eta3)/8;
    gradPhi_Eta[1].z() = -(1 + eta1)*(1 - eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex C (cellVertices[2])
    // ---------------------------------------------------------------------
    gradPhi_Eta[2].x() = -(1 + eta2)*(1 - eta3)/8;
    gradPhi_Eta[2].y() =  (1 - eta1)*(1 - eta3)/8;
    gradPhi_Eta[2].z() = -(1 - eta1)*(1 + eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex D (cellVertices[3])
    // ---------------------------------------------------------------------
    gradPhi_Eta[3].x() =  (1 + eta2)*(1 - eta3)/8;
    gradPhi_Eta[3].y() =  (1 + eta1)*(1 - eta3)/8;
    gradPhi_Eta[3].z() = -(1 + eta1)*(1 + eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex E (cellVertices[4])
    // ---------------------------------------------------------------------
    gradPhi_Eta[4].x() = 0;
    gradPhi_Eta[4].y() = 0;
    gradPhi_Eta[4].z() = 0.5;

    // ---------------------------------------------------------------------
    // Compute physical gradients: d(x,y,z)/d(eta1,eta2,eta3)
    // ---------------------------------------------------------------------
    vector gradX1_Eta(Zero); // d(x)/d(eta)
    vector gradX2_Eta(Zero); // d(y)/d(eta)
    vector gradX3_Eta(Zero); // d(z)/d(eta)

    for (label i = 0; i < cellVertices.size(); ++i)
    {
        gradX1_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].x();
        gradX1_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].x();
        gradX1_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].x();

        gradX2_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].y();
        gradX2_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].y();
        gradX2_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].y();

        gradX3_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].z();
        gradX3_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].z();
        gradX3_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].z();
    }

    // ---------------------------------------------------------------------
    // Jacobian tensor (columns = gradients of x, y, z w.r.t eta)
    // ---------------------------------------------------------------------
    const tensor J3D
    (
        gradX1_Eta,
        gradX2_Eta,
        gradX3_Eta
    );

    // ---------------------------------------------------------------------
    // Return determinant of Jacobian
    // ---------------------------------------------------------------------
    return det(J3D);
}

Foam::scalar Foam::TetraJacobian3D
(
    const vector& gaussPt,             // Gauss point in reference space (eta)
    const List<vector>& cellVertices   // Physical coordinates of tetrahedron vertices
)
{
    // ---------------------------------------------------------------------
    // Extract reference coordinates (eta-space)
    // ---------------------------------------------------------------------
    const scalar eta1 = gaussPt.x();
    const scalar eta2 = gaussPt.y();
    const scalar eta3 = gaussPt.z();

    // ---------------------------------------------------------------------
    // Gradients of shape functions w.r.t eta (dPhi/dEta)
    //
    // Vertex correspondence:
    //
    //   cellVertices id -> tetrahedron vertex -> gradPhi_Eta id
    //   ---------------------------------------------------
    //        0        ->        A         ->        0
    //        1        ->        B         ->        1
    //        2        ->        C         ->        2
    //        3        ->        D         ->        3
    //
    // ---------------------------------------------------------------------
    List<vector> gradPhi_Eta(cellVertices.size(), Zero);

    // ---------------------------------------------------------------------
    // Vertex A (cellVertices[0])
    // ---------------------------------------------------------------------
    gradPhi_Eta[0].x() = -(1 - eta2)*(1 - eta3)/8;
    gradPhi_Eta[0].y() = -(1 - eta1)*(1 - eta3)/8;
    gradPhi_Eta[0].z() = -(1 - eta1)*(1 - eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex B (cellVertices[1])
    // ---------------------------------------------------------------------
    gradPhi_Eta[1].x() =  (1 - eta2)*(1 - eta3)/8;
    gradPhi_Eta[1].y() = -(1 + eta1)*(1 - eta3)/8;
    gradPhi_Eta[1].z() = -(1 + eta1)*(1 - eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex C (cellVertices[2])
    // ---------------------------------------------------------------------
    gradPhi_Eta[2].x() = -(1 + eta2)*(1 - eta3)/8;
    gradPhi_Eta[2].y() =  (1 - eta1)*(1 - eta3)/8;
    gradPhi_Eta[2].z() = -(1 - eta1)*(1 + eta2)/8;

    // ---------------------------------------------------------------------
    // Vertex D (cellVertices[3])
    // ---------------------------------------------------------------------
    gradPhi_Eta[3].x() = 0;
    gradPhi_Eta[3].y() = 0;
    gradPhi_Eta[3].z() = 0.5;

    // ---------------------------------------------------------------------
    // Compute physical gradients: d(x,y,z)/d(eta1,eta2,eta3)
    // ---------------------------------------------------------------------
    vector gradX1_Eta(Zero); // d(x)/d(eta)
    vector gradX2_Eta(Zero); // d(y)/d(eta)
    vector gradX3_Eta(Zero); // d(z)/d(eta)

    for (label i = 0; i < cellVertices.size(); ++i)
    {
        gradX1_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].x();
        gradX1_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].x();
        gradX1_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].x();

        gradX2_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].y();
        gradX2_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].y();
        gradX2_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].y();

        gradX3_Eta.x() += gradPhi_Eta[i].x()*cellVertices[i].z();
        gradX3_Eta.y() += gradPhi_Eta[i].y()*cellVertices[i].z();
        gradX3_Eta.z() += gradPhi_Eta[i].z()*cellVertices[i].z();
    }

    // ---------------------------------------------------------------------
    // Jacobian tensor (columns = gradients of x, y, z w.r.t eta)
    // ---------------------------------------------------------------------
    const tensor J3D
    (
        gradX1_Eta,
        gradX2_Eta,
        gradX3_Eta
    );

    // ---------------------------------------------------------------------
    // Return determinant of Jacobian
    // ---------------------------------------------------------------------
    return det(J3D);
}

} // namespace math
} // namespace Foam
// ************************************************************************* //
