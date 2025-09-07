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

#include "dgMath.H"
#include "error.H"

namespace Foam
{
namespace math
{

// Construct scalar integral from Gauss samples
Foam::scalar gaussIntegral
(
    const List<scalar>& values,
    const List<scalar>& weights
)
{
    // Validate sizes
    if (values.size() != weights.size())
    {
        FatalErrorInFunction
            << "Size mismatch: values.size() = " << values.size()
            << ", weights.size() = " << weights.size() << nl
            << exit(FatalError);
    }

    // Accumulate sum_i values[i] * weights[i]
    scalar I = 0.0;
    forAll(values, i)
    {
        I += values[i] * weights[i];
    }
    return I;
}


// Construct vector integral from Gauss samples (component-wise)
Foam::vector gaussIntegral
(
    const List<vector>& values,
    const List<scalar>& weights
)
{
    // Validate sizes
    if (values.size() != weights.size())
    {
        FatalErrorInFunction
            << "Size mismatch: values.size() = " << values.size()
            << ", weights.size() = " << weights.size() << nl
            << exit(FatalError);
    }

    // Accumulate sum_i values[i] * weights[i]
    vector I(Zero);
    forAll(values, i)
    {
        I += values[i] * weights[i];  // vector * scalar
    }
    return I;
}


// Construct tensor integral from Gauss samples (component-wise)
Foam::tensor gaussIntegral
(
    const List<tensor>& values,
    const List<scalar>& weights
)
{
    // Validate sizes
    if (values.size() != weights.size())
    {
        FatalErrorInFunction
            << "Size mismatch: values.size() = " << values.size()
            << ", weights.size() = " << weights.size() << nl
            << exit(FatalError);
    }

    // Accumulate sum_i values[i] * weights[i]
    tensor I(Zero);
    forAll(values, i)
    {
        I += values[i] * weights[i];  // tensor * scalar
    }
    return I;
}

} // namespace math
} // namespace Foam
// ************************************************************************* //
