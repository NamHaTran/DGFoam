/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
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

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgCompressibleBoundaryState_H
#define Foam_dgCompressibleBoundaryState_H

#include "label.H"
#include "scalar.H"
#include "tensor.H"
#include "vector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Struct dgCompressibleBoundaryState Declaration
\*---------------------------------------------------------------------------*/

/**
 * \struct Foam::dgCompressibleBoundaryState
 * \brief Conservative and differential state at one boundary Gauss point.
 *
 * \details
 * The structure gathers all quantities needed by reusable compressible
 * velocity-slip and temperature-jump models.  Values with the \c Minus suffix
 * are reconstructed from the owner cell side of the boundary, while values
 * with the \c Plus suffix are the current ghost/exterior state.
 */
struct dgCompressibleBoundaryState
{
    /// Owner cell label.
    label cellID;

    /// Local face label within the owner cell.
    label faceLocalID;

    /// Local Gauss point label on the face.
    label localGauss;

    /// Unit face normal used by the boundary flux.
    vector n;

    /// Normal distance coefficient between interior and boundary state.
    scalar normalDeltaCoeff;

    /// Interior density.
    scalar rhoMinus;

    /// Interior momentum.
    vector rhoUMinus;

    /// Interior total energy per volume.
    scalar EMinus;

    /// Exterior/ghost density.
    scalar rhoPlus;

    /// Exterior/ghost momentum.
    vector rhoUPlus;

    /// Exterior/ghost total energy per volume.
    scalar EPlus;

    /// Interior velocity-gradient tensor.
    tensor gradUMinus;

    /// Interior temperature-gradient vector.
    vector gradTMinus;

    /// Construct an empty, invalid-index boundary state.
    dgCompressibleBoundaryState()
    :
        cellID(-1),
        faceLocalID(-1),
        localGauss(-1),
        n(Zero),
        normalDeltaCoeff(Zero),
        rhoMinus(Zero),
        rhoUMinus(Zero),
        EMinus(Zero),
        rhoPlus(Zero),
        rhoUPlus(Zero),
        EPlus(Zero),
        gradUMinus(Zero),
        gradTMinus(Zero)
    {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
