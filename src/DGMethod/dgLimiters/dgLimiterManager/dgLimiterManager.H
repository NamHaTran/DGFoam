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
    Foam::dgLimiterManager

Description
    Manager for DG limiters configured in constant/dgOptions.

SourceFiles
    dgLimiterManager.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgLimiterManager_H
#define Foam_dgLimiterManager_H

#include "IOdictionary.H"
#include "PtrList.H"
#include "List.H"
#include "word.H"
#include "dgGeomMesh.H"
#include "dgLimiter.H"

namespace Foam
{

class Ostream;

/*---------------------------------------------------------------------------*\
                     Class dgLimiterManager Declaration
\*---------------------------------------------------------------------------*/

class dgLimiterManager
:
    public IOdictionary
{
private:

        /** Shared DG geometric mesh. */
        const dgGeomMesh& mesh_;

        /** Runtime-selected limiter instances. */
        PtrList<dgLimiter> limiters_;

        /** Limiter names following the ordering in \c limiters_. */
        List<word> limiterNames_;

        /** Create the IOobject for \c constant/dgOptions. */
        IOobject createIOobject(const dgGeomMesh& mesh) const;

        /** Construct the limiter list from the current dictionary content. */
        void constructLimiters();

        //- No copy construct
        dgLimiterManager(const dgLimiterManager&) = delete;

        //- No copy assignment
        void operator=(const dgLimiterManager&) = delete;

public:

        /**
         * \brief Construct and read limiter configuration from \c constant/dgOptions.
         *
         * If the dictionary is absent, the manager remains empty and
         * \c correct() becomes a no-op.
         */
        explicit dgLimiterManager(const dgGeomMesh& mesh);

        //- Destructor
        ~dgLimiterManager() = default;

        /** Return the number of configured limiters. */
        inline label size() const
        {
            return limiters_.size();
        }

        /** Return true when no limiter is configured. */
        inline bool empty() const
        {
            return limiters_.empty();
        }

        /** Return the configured limiter names. */
        inline const List<word>& names() const
        {
            return limiterNames_;
        }

        /** Return true if a limiter with the given name exists. */
        bool has(const word& name) const;

        /** Return const access to a limiter by name. */
        const dgLimiter& limiter(const word& name) const;

        /** Return mutable access to a limiter by name. */
        dgLimiter& limiter(const word& name);

        /** Apply every configured limiter in dictionary order. */
        void correct();

        /** Re-read \c dgOptions and rebuild the limiter list if modified. */
        virtual bool read() override;

        /** Print the configured limiter names and limiter types. */
        void listLimiters(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
