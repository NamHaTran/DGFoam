/**
 * \brief Point-wise conservative and derived quantities written to VTK.
 *
 * \details
 * The export utility reconstructs DG modal fields at every VTK Lagrange node.
 * For each reconstructed point it stores:
 * - conservative variables \f$\rho\f$, \f$\rho \mathbf{U}\f$, \f$E\f$,
 * - primitive variables \f$\mathbf{U}\f$, \f$p\f$, \f$T\f$,
 * - and auxiliary thermodynamic data \f$e\f$, \f$a\f$.
 *
 * The arrays are accumulated in the same order as `grid.points`, so they can
 * be passed directly to `dgVtkLagrange::writeVtuPiece()`.
 */
struct ExportPointData
{
    std::vector<scalar> rho;
    std::vector<vector> rhoU;
    std::vector<scalar> E;
    std::vector<vector> U;
    std::vector<scalar> e;
    std::vector<scalar> p;
    std::vector<scalar> T;
    std::vector<scalar> a;
};


/** \brief Conservative state sampled at one Lagrange node. */
struct ConservativeSample
{
    scalar rho;
    vector rhoU;
    scalar E;
};


constexpr unsigned char vtkLagrangeTriangleType = 69u;
constexpr unsigned char vtkLagrangeQuadType = 70u;


inline label trianglePointCount(const label order)
{
    return (order + 1)*(order + 2)/2;
}


inline label trianglePointIndex
(
    const label order,
    const label row,
    const label col
)
{
    return row*(order + 1) - row*(row - 1)/2 + col;
}


inline vector referenceVertexEta
(
    const dgCellType type,
    const label vertexI
)
{
    switch (type)
    {
        case dgCellType::HEX:
        {
            static const vector eta[8] =
            {
                vector(-1, -1, -1),
                vector( 1, -1, -1),
                vector( 1,  1, -1),
                vector(-1,  1, -1),
                vector(-1, -1,  1),
                vector( 1, -1,  1),
                vector( 1,  1,  1),
                vector(-1,  1,  1)
            };

            return eta[vertexI];
        }

        case dgCellType::PRISM:
        {
            static const vector eta[6] =
            {
                vector(-1, -1,  1),
                vector(-1, -1, -1),
                vector( 1, -1, -1),
                vector(-1,  1,  1),
                vector(-1,  1, -1),
                vector( 1,  1, -1)
            };

            return eta[vertexI];
        }

        case dgCellType::PYRAMID:
        {
            static const vector eta[5] =
            {
                vector(-1, -1, -1),
                vector( 1, -1, -1),
                vector( 1,  1, -1),
                vector(-1,  1, -1),
                vector(-1, -1,  1)
            };

            return eta[vertexI];
        }

        case dgCellType::TET:
        {
            static const vector eta[4] =
            {
                vector(-1, -1, -1),
                vector( 1, -1, -1),
                vector(-1,  1, -1),
                vector(-1, -1,  1)
            };

            return eta[vertexI];
        }

        default:
            FatalErrorInFunction
                << "Unsupported DG cell type " << type
                << exit(FatalError);
    }

    return vector::zero;
}


inline label findOwnerVertex
(
    const point& faceVertex,
    const pointField& ownerCellVertices
)
{
    label nearestVertexI = -1;
    scalar nearestDist = GREAT;

    forAll(ownerCellVertices, vertexI)
    {
        const scalar dist = mag(faceVertex - ownerCellVertices[vertexI]);

        if (dist < nearestDist)
        {
            nearestDist = dist;
            nearestVertexI = vertexI;
        }
    }

    const scalar tol = 1e-8*max(scalar(1), mag(faceVertex));

    if (nearestVertexI < 0 || nearestDist > tol)
    {
        FatalErrorInFunction
            << "Cannot match patch vertex " << faceVertex
            << " to an owner-cell vertex. Nearest distance = "
            << nearestDist << ", tolerance = " << tol
            << exit(FatalError);
    }

    return nearestVertexI;
}


inline vector referencePointEta
(
    const dgCellType type,
    const UList<label>& ownerVertexIDs,
    const UList<scalar>& weights
)
{
    if (type == dgCellType::TET)
    {
        scalar lambda[4] = {0, 0, 0, 0};

        forAll(ownerVertexIDs, i)
        {
            lambda[ownerVertexIDs[i]] += weights[i];
        }

        const scalar lamB = lambda[1];
        const scalar lamC = lambda[2];
        const scalar lamD = lambda[3];

        const scalar eta3 = 2.0*lamD - 1.0;
        const scalar denom2 = 1.0 - lamD;

        if (mag(denom2) <= SMALL)
        {
            return vector(-1.0, -1.0, eta3);
        }

        const scalar eta2 = 2.0*lamC/denom2 - 1.0;
        const scalar denom1 = 1.0 - lamC - lamD;

        if (mag(denom1) <= SMALL)
        {
            return vector(-1.0, eta2, eta3);
        }

        const scalar eta1 = 2.0*lamB/denom1 - 1.0;

        return vector(eta1, eta2, eta3);
    }

    vector eta = vector::zero;

    forAll(ownerVertexIDs, i)
    {
        eta += weights[i]*referenceVertexEta(type, ownerVertexIDs[i]);
    }

    return eta;
}


inline label vtkTriangleIndex
(
    const std::array<label, 3>& bindex,
    label order
)
{
    label index = 0;
    label maxCoord = order;
    label minCoord = 0;
    const label bmin = min(min(bindex[0], bindex[1]), bindex[2]);

    while (bmin > minCoord)
    {
        index += 3*order;
        maxCoord -= 2;
        ++minCoord;
        order -= 3;
    }

    for (label dim = 0; dim < 3; ++dim)
    {
        if (bindex[(dim + 2)%3] == maxCoord)
        {
            return index;
        }

        ++index;
    }

    for (label dim = 0; dim < 3; ++dim)
    {
        if (bindex[(dim + 1)%3] == minCoord)
        {
            return index + bindex[dim] - (minCoord + 1);
        }

        index += maxCoord - (minCoord + 1);
    }

    return index;
}


std::vector<label> vtkLagrangeTriangleNodemap(const label order)
{
    std::vector<label> nodemap(trianglePointCount(order), -1);

    for (label row = 0; row <= order; ++row)
    {
        for (label col = 0; col <= order - row; ++col)
        {
            const label l0 = order - row - col;
            const label l1 = col;
            const label l2 = row;
            const std::array<label, 3> bindex = {l1, l2, l0};
            const label vtkPointI = vtkTriangleIndex(bindex, order);

            nodemap[vtkPointI] = trianglePointIndex(order, row, col);
        }
    }

    for (const label pointI : nodemap)
    {
        if (pointI < 0)
        {
            FatalErrorInFunction
                << "Failed to construct a complete VTK triangle nodemap for "
                << "order " << order << exit(FatalError);
        }
    }

    return nodemap;
}


std::vector<label> vtkLagrangeQuadNodemap(const label order)
{
    std::vector<label> nodemap;
    nodemap.reserve((order + 1)*(order + 1));

    nodemap.push_back(0);
    nodemap.push_back(order);
    nodemap.push_back(order*(order + 1) + order);
    nodemap.push_back(order*(order + 1));

    for (label col = 1; col < order; ++col)
    {
        nodemap.push_back(col);
    }

    for (label row = 1; row < order; ++row)
    {
        nodemap.push_back(row*(order + 1) + order);
    }

    for (label col = order - 1; col > 0; --col)
    {
        nodemap.push_back(order*(order + 1) + col);
    }

    for (label row = order - 1; row > 0; --row)
    {
        nodemap.push_back(row*(order + 1));
    }

    for (label row = 1; row < order; ++row)
    {
        for (label col = 1; col < order; ++col)
        {
            nodemap.push_back(row*(order + 1) + col);
        }
    }

    return nodemap;
}


inline void computeOwnerBasisAt
(
    const dgCellType type,
    const vector& eta,
    const label pOrder,
    List<scalar>& basis
)
{
    List<scalar> dEta1;
    List<scalar> dEta2;
    List<scalar> dEta3;

    switch (type)
    {
        case dgCellType::HEX:
            computeHexBasisAndDerivatives
            (
                eta,
                pOrder,
                basis,
                dEta1,
                dEta2,
                dEta3
            );
            break;

        case dgCellType::PRISM:
            computePrismBasisAndDerivatives
            (
                eta,
                pOrder,
                basis,
                dEta1,
                dEta2,
                dEta3
            );
            break;

        case dgCellType::PYRAMID:
            computePyramidBasisAndDerivatives
            (
                eta,
                pOrder,
                basis,
                dEta1,
                dEta2,
                dEta3
            );
            break;

        case dgCellType::TET:
            computeTetBasisAndDerivatives
            (
                eta,
                pOrder,
                basis,
                dEta1,
                dEta2,
                dEta3
            );
            break;

        default:
            FatalErrorInFunction
                << "Unsupported DG cell type " << type
                << exit(FatalError);
    }
}
