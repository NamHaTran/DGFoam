/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2026 Ha Nam Tran
-------------------------------------------------------------------------------
License
    This file is part of DGFoam.

    DGFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::dgPDETerm

Description
    Mesh-level polymorphic interface for DG PDE terms.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dgPDETerm_H
#define Foam_dgPDETerm_H

#include "GaussField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dgPDETerm Declaration
\*---------------------------------------------------------------------------*/

/**
 * \class Foam::dgPDETerm
 * \brief Common interface for DG flux and source terms assembled by
 *        dgResidual.
 *
 * \tparam UType Type of the conservative variable and residual entry.
 *               Current use cases are \c scalar for mass/energy and
 *               \c vector for momentum.
 * \tparam FType Type of the divergence flux associated with \c UType.
 *               Current use cases are \c vector for scalar equations and
 *               \c tensor for vector equations.
 *
 * \details
 * A \c dgPDETerm does not assemble a residual by itself.  It exposes either
 * a flux field, a source field, or both, and \c dgResidual performs the common
     * integration for all terms belonging to one equation on one active cell.
 *
     * The interface is intentionally reusable.  A concrete term is constructed
     * once from the mesh, receives its active cell fields through \c add(), and
     * can then be passed to \c dgResidual for that cell.
 *
 * Flux terms such as \c dgConvectiveTerm and \c dgDiffusiveTerm override
 * \c fluxFieldPtr() and usually \c updateFlux().  Source terms such as
 * \c dgSourceTerm override \c sourceFieldPtr() and leave \c updateFlux() as a
 * no-op.
 */
template<class UType, class FType>
class dgPDETerm
{
public:

    /** Type of the equation unknown and residual entry. */
    using solution_type = UType;

    /** Type of the divergence flux integrated against \f$\nabla\phi\f$. */
    using flux_type = FType;

    /** Virtual destructor for safe use through base references/pointers. */
    virtual ~dgPDETerm() = default;

    /**
     * \brief Return true when this term represents a convective flux.
     *
     * This flag is diagnostic/classification metadata.  Residual assembly is
     * driven by whether \c fluxFieldPtr() or \c sourceFieldPtr() returns a
     * non-null pointer.
     */
    virtual bool isConvection() const { return false; }

    /**
     * \brief Return true when this term represents a diffusive flux.
     *
     * This includes physical diffusion and artificial diffusion terms.
     */
    virtual bool isDiffusion() const { return false; }

    /**
     * \brief Return true when this term represents a volume source.
     */
    virtual bool isSource() const { return false; }

    /**
     * \brief Update numerical face fluxes for flux-like terms on one cell.
     *
     * The base implementation is a no-op so source-only terms can be stored in
     * the same list without special handling.
     */
    virtual void updateFlux(const label cellID) {}

    /**
     * \brief Return the term flux field, or \c nullptr for source-only terms.
     *
     * The returned field contains cell-interior flux values and face numerical
     * flux values.  \c dgResidual uses the cell field for the volume term and
     * the face field for the surface term.
     */
    virtual const GaussField<FType>* fluxFieldPtr() const { return nullptr; }

    /**
     * \brief Return the term source field, or \c nullptr for flux-only terms.
     *
     * Source fields are integrated as \f$\int_K \phi_i S\,dV\f$ and are not
     * treated as divergence fluxes.
     */
    virtual const cellGaussField<UType>* sourceFieldPtr() const { return nullptr; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
