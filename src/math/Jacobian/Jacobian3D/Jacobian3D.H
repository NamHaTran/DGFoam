/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Description
    A collection of utilities for evaluating 3D Jacobian

See also

SourceFiles
    Jacobian3D.C

\*---------------------------------------------------------------------------*/

#ifndef Jacobian3D_H
#define Jacobian3D_H

#include "dgCellType.H"
#include "scalar.H"
#include "Tensor.H"
#include "label.H"
#include "List.H"
#include "error.H"  // For FatalErrorInFunction

namespace Foam
{
namespace math
{

Foam::scalar Foam::HexJacobian3D
(
    const vector& gaussPt,             // Gauss point in reference space (eta)
    const List<vector>& cellVertices   // Physical coordinates of hex vertices
);

Foam::scalar Foam::PrismJacobian3D
(
    const vector& gaussPt,             // Gauss point in reference space (eta)
    const List<vector>& cellVertices   // Physical coordinates of prism vertices
);

Foam::scalar Foam::PyramidJacobian3D
(
    const vector& gaussPt,             // Gauss point in reference space (eta)
    const List<vector>& cellVertices   // Physical coordinates of pyramid vertices
);

Foam::scalar Foam::TetraJacobian3D
(
    const vector& gaussPt,             // Gauss point in reference space (eta)
    const List<vector>& cellVertices   // Physical coordinates of tetrahedron vertices
);

} // namespace math
} // namespace Foam

#endif

// ************************************************************************* //

