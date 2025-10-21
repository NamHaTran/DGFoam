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

Class
    Foam::dgFluxSolverHelper

Description
    dgFluxSolverHelper is a per-cell cache for intermediate states
    (e.g. Roe, HLLC, HLLE) used in flux solvers. It stores one state
    per Gauss point across all faces of a cell, and allows flux solvers
    to reuse previously computed data to avoid redundant calculations.

    The class is designed to be easily extensible when adding support
    for new flux solvers. To add a new intermediate state (e.g. Roe), follow
    these steps:

    1) Add a new List container to store the state per Gauss point:
        Example:
            List<RoeState> RoeStateList_;

    2) Add new methods to interface with the state:
        - setRoeAvailable()          // To mark that Roe solver is in use
        - accessRoeState(label id)   // To get mutable access to RoeState
        - accessRoeState(label id) const // To get const access

    3) Add a method to reset the validity of the new state:
        - resetRoeStateValidity()

       Optionally, add this reset step into the global reset() method.

    Each intermediate state struct must include a `bool valid` flag
    to indicate whether the data has already been computed for the
    current time step.

SourceFiles
    dgFluxSolverHelper.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgFluxSolverHelper_H
#define Foam_dgFluxSolverHelper_H

#include "HLLCState.H"
#include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dgFluxSolverHelper Declaration
\*---------------------------------------------------------------------------*/

class dgFluxSolverHelper
{
private:

    //- Number of Gauss points across all faces in the cell
    label nGauss_;

    //- Flag: is HLLC solver used for this cell?
    bool hasHLLC_ = false;

    //- HLLC intermediate state per Gauss point
    List<HLLCState> HLLCStateList_;

public:

    //- Constructor with number of Gauss points
    dgFluxSolverHelper(const label nGauss)
    :
        nGauss_(nGauss),
        HLLCStateList_(nGauss)
    {}

    //- Set that HLLC solver is active
    inline void setHLLCAvailable()
    {
        hasHLLC_ = true;
    }

    //- Return non-const pointer to HLLCState at given gaussID
    inline HLLCState* accessHLLCState(const label gaussID)
    {
        return &HLLCStateList_[gaussID];
    }

    //- Return const pointer to HLLCState at given gaussID
    inline const HLLCState* accessHLLCState(const label gaussID) const
    {
        return &HLLCStateList_[gaussID];
    }

    //- Reset valid flag for all HLLC states
    inline void resetHLLCStateValidity()
    {
        forAll(HLLCStateList_, i)
        {
            HLLCStateList_[i].valid = false;
        }
    }

    //- Global reset (for future extension)
    inline void reset()
    {
        resetHLLCStateValidity();
        hasHLLC_ = false;
    }

    //- Return number of Gauss points
    inline label nGauss() const
    {
        return nGauss_;
    }
};

} // namespace Foam

#endif

// ************************************************************************* //