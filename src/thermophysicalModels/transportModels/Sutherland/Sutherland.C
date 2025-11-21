/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright
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

#include "Sutherland.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Sutherland, 0);
addToRunTimeSelectionTable(transportLaw, Sutherland, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Sutherland::Sutherland
(
    const word& name,
    const dictionary& dict,
    const dgGeomMesh& mesh,
    const thermoLaw& thermo
)
:
    transportLaw(name, dict, mesh, thermo),
    As_(1.458e-5),
    S_(110.4),
    Pr0_(0.72)
{
    read();
}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Sutherland::calcMu
(
    const label cellI,
    const GaussField<scalar>& T,
    GaussField<scalar>& mu
) const
{
    // μ = As * T^(3/2) / (T + S)

    // Guard: T must be positive
    if (T <= scalar(0))
    {
        FatalErrorInFunction
            << "Non-positive temperature in Sutherland::calcMu()" << nl
            << exit(FatalError);
    }

    GaussField<scalar> denom(cellI, &mesh_);
    denom = T + S_;   // S_ is scalar, auto-broadcast

    GaussField<scalar> T32(cellI, &mesh_);
    T32 = pow(T, scalar(1.5));

    mu = (As_ * T32) / denom;
}

void Sutherland::calcKappa
(
    const label cellI,
    const GaussField<scalar>& T,
    GaussField<scalar>& kappa
) const
{
    // classical: k = μ Cp / Pr
    GaussField<scalar> mu(cellI, &mesh_);
    calcMu(cellI, T, mu);

    GaussField<scalar> Cp(cellI, &mesh_);
    thermo_.calcCp(cellI, T, Cp);

    // Pr constant
    const scalar Pr = Pr0_;

    kappa = (mu * Cp) / Pr;
}

void Sutherland::calcPr
(
    const label cellI,
    const GaussField<scalar>& T,
    GaussField<scalar>& Pr
) const
{
    Pr = Pr0_;   // broadcast assignment to all Gauss points
}

void Sutherland::read()
{
    if (coeff_.found("As")) { As_  = readScalar(coeff_.lookup("As")); }
    if (coeff_.found("S"))  { S_   = readScalar(coeff_.lookup("S")); }
    if (coeff_.found("Pr")) { Pr0_ = readScalar(coeff_.lookup("Pr")); }

    if (As_ <= 0)
    {
        FatalErrorInFunction
            << "Invalid Sutherland 'As'. Must be > 0." << exit(FatalError);
    }
    if (S_ <= 0)
    {
        FatalErrorInFunction
            << "Invalid Sutherland 'S'. Must be > 0." << exit(FatalError);
    }
    if (Pr0_ <= 0)
    {
        FatalErrorInFunction
            << "Invalid Prandtl number. Must be > 0." << exit(FatalError);
    }
}

} // End namespace Foam
