/**
 * \brief Reconstruct a modal DG field at one interpolation point in one cell.
 *
 * \tparam Type Scalar or vector field value type.
 * \param field DG field containing modal coefficients on every cell.
 * \param cellI Cell index whose local modal expansion is evaluated.
 * \param basis Values of the local basis functions at the target point.
 * \return Reconstructed physical value \f$\sum_i \phi_i(\xi)\hat{u}_i\f$.
 *
 * \details
 * The DG solution is stored as modal coefficients per cell. This helper
 * evaluates the polynomial inside `cellI` by forming the dot product between:
 * - the modal coefficients in `field.dof()[cellI]`, and
 * - the basis values already evaluated at the local export point.
 *
 * A strict size check is performed before reconstruction because a mismatch
 * would indicate inconsistent polynomial order or an invalid nodal basis.
 */
template<class Type>
Type evaluateModalField
(
    const dgField<Type>& field,
    const label cellI,
    const List<scalar>& basis
)
{
    const cellDof<Type>& cellModes = field.dof()[cellI];

    if (cellModes.nDof() != basis.size())
    {
        FatalErrorInFunction
            << "Basis size mismatch while evaluating field " << field.name()
            << " in cell " << cellI
            << ". DoFs = " << cellModes.nDof()
            << ", basis size = " << basis.size()
            << exit(FatalError);
    }

    Type value = pTraits<Type>::zero;

    for (label modeI = 0; modeI < cellModes.nDof(); ++modeI)
    {
        value += basis[modeI]*cellModes[modeI];
    }

    return value;
}


/**
 * \brief Return the zeroth modal coefficient for one DG field in one cell.
 *
 * \details
 * The modal basis used by DGFoam stores the cell-average value in DoF0. This
 * value is used as a constant reconstruction when theta-based smoothing clamps
 * a cell to its mean value.
 */
template<class Type>
Type evaluateCellMeanField
(
    const dgField<Type>& field,
    const label cellI
)
{
    const cellDof<Type>& cellModes = field.dof()[cellI];

    if (cellModes.nDof() == 0)
    {
        FatalErrorInFunction
            << "Field " << field.name() << " has no DoF in cell "
            << cellI << exit(FatalError);
    }

    return cellModes[0];
}


/**
 * \brief Zhang-Shu limiter coefficient clamped to [0, 1].
 */
scalar limiterTheta
(
    const scalar meanValue,
    const scalar minValue,
    const scalar omega,
    const scalar epsilon
)
{
    const scalar scale = max(scalar(1), mag(meanValue));

    if (mag(meanValue - minValue) <= epsilon*scale)
    {
        return 1.0;
    }

    const scalar rawTheta = (meanValue - omega)/(meanValue - minValue);
    return min(max(rawTheta, scalar(0)), scalar(1));
}


/**
 * \brief Return rho*e from conservative variables with a guarded rho.
 */
scalar thermoEnergyDensity
(
    const scalar rho,
    const vector& rhoU,
    const scalar E,
    const scalar epsilon
)
{
    const scalar rhoSafe =
        (mag(rho) > epsilon ? rho : (rho >= 0 ? epsilon : -epsilon));

    return E - 0.5*magSqr(rhoU)/rhoSafe;
}


/**
 * \brief Scale the non-mean modal content in one cell.
 */
template<class Type>
void scaleHighOrderModes
(
    dgField<Type>& field,
    const label cellI,
    const scalar theta
)
{
    if (!field.hasDof() || theta >= 1.0)
    {
        return;
    }

    cellDof<Type>& cellModes = field.dof()[cellI];

    for (label modeI = 1; modeI < cellModes.nDof(); ++modeI)
    {
        cellModes[modeI] *= theta;
    }

    field.dof().updateCellDof(cellI);
}


/**
 * \brief Clamp conservative fields to DoF0 using stored limiter theta.
 *
 * Cells satisfying the theta condition and their local face-neighbour cells are
 * clamped.
 */
void clampConservativeFieldsByTheta
(
    const dgGeomMesh& dgMesh,
    dgField<scalar>& rhoField,
    dgField<vector>& rhoUField,
    dgField<scalar>& EField,
    const scalar clampTheta,
    const volScalarField& theta1Field,
    const volScalarField& theta2Field
)
{
    List<bool> clampCells(dgMesh.nCells(), false);

    for (label cellI = 0; cellI < dgMesh.nCells(); ++cellI)
    {
        const dgGeomCell* cellPtr = dgMesh.cells()[cellI];

        if (!cellPtr)
        {
            continue;
        }

        const scalar theta1 = theta1Field.primitiveField()[cellI];
        const scalar theta2 = theta2Field.primitiveField()[cellI];
        const bool clampCell =
            (theta1 >= clampTheta && theta1 < scalar(1))
         || (theta2 >= clampTheta && theta2 < scalar(1));

        if (clampCell)
        {
            clampCells[cellI] = true;

            const labelList& neighbourCells = cellPtr->neighborCells();

            forAll(neighbourCells, neighbourI)
            {
                const label neighbourCellI = neighbourCells[neighbourI];

                if
                (
                    neighbourCellI >= 0
                 && neighbourCellI < dgMesh.nCells()
                 && dgMesh.cells()[neighbourCellI]
                )
                {
                    clampCells[neighbourCellI] = true;
                }
            }
        }
    }

    label nClampedCells = 0;

    forAll(clampCells, cellI)
    {
        if (clampCells[cellI])
        {
            scaleHighOrderModes(rhoField, cellI, scalar(0));
            scaleHighOrderModes(rhoUField, cellI, scalar(0));
            scaleHighOrderModes(EField, cellI, scalar(0));
            ++nClampedCells;
        }
    }

    reduce(nClampedCells, sumOp<label>());

    Info<< "smoothByClampingTheta " << clampTheta
        << ": clamped rho/rhoU/E high-order modes to zero in marked cells "
        << "and their local neighbours, total "
        << nClampedCells << " cell(s)." << nl;
}


/**
 * \brief Limit modal conservative fields using the same Lagrange nodes that are
 *        later written to VTK.
 */
void limitConservativeFieldsAtLagrangeNodes
(
    const dgGeomMesh& dgMesh,
    dgField<scalar>& rhoField,
    dgField<vector>& rhoUField,
    dgField<scalar>& EField,
    const scalar epsilon
)
{
    label nLimitedCells = 0;

    const label pOrder = dgMesh.pOrder();
    const label hoOrder = dgVtkLagrange::exportOrder(pOrder);

    if (hoOrder > 8)
    {
        FatalErrorInFunction
            << "Current dgToVTKLagrange implementation supports "
            << "high-order export up to order 8. "
            << "Detected DG pOrder = " << pOrder
            << ", which maps to VTK order " << hoOrder
            << exit(FatalError);
    }

    for (label cellI = 0; cellI < dgMesh.nCells(); ++cellI)
    {
        const dgGeomCell* cellPtr = dgMesh.cells()[cellI];

        if (!cellPtr)
        {
            continue;
        }

        if (rhoField.dof()[cellI].nDof() <= 1)
        {
            continue;
        }

        const dgCellType type = cellPtr->type();

        if (!dgVtkLagrange::supportsNativeLagrange(type))
        {
            if (type == dgCellType::PYRAMID)
            {
                continue;
            }

            FatalErrorInFunction
                << "Unsupported DG cell type "
                << dgVtkLagrange::cellTypeName(type)
                << " in cell " << cellI << exit(FatalError);
        }

        const std::vector<vector> pyfrPts =
            dgVtkLagrange::pyfrStdElementPoints(type, hoOrder);
        const std::vector<label> nodemap =
            dgVtkLagrange::vtkNodemap(type, pyfrPts.size());

        DynamicList<ConservativeSample> samples(nodemap.size());
        scalar minRho = GREAT;

        for (const label pyfrPointI : nodemap)
        {
            const vector eta =
                dgVtkLagrange::xiToEta(type, pyfrPts[pyfrPointI]);

            List<scalar> basis;
            dgVtkLagrange::computeBasisAt(type, eta, pOrder, basis);

            ConservativeSample sample;
            sample.rho = evaluateModalField(rhoField, cellI, basis);
            sample.rhoU = evaluateModalField(rhoUField, cellI, basis);
            sample.E = evaluateModalField(EField, cellI, basis);

            minRho = min(minRho, sample.rho);
            samples.append(sample);
        }

        const scalar meanRho = rhoField.dof()[cellI].dof()[0];
        const vector& meanRhoU = rhoUField.dof()[cellI].dof()[0];
        const scalar meanE = EField.dof()[cellI].dof()[0];

        const scalar theta1 =
            limiterTheta(meanRho, minRho, min(epsilon, meanRho), epsilon);

        scalar minRhoThermoEnergy = GREAT;

        forAll(samples, sampleI)
        {
            const ConservativeSample& sample = samples[sampleI];
            const scalar limitedRho =
                meanRho + theta1*(sample.rho - meanRho);

            minRhoThermoEnergy =
                min
                (
                    minRhoThermoEnergy,
                    thermoEnergyDensity
                    (
                        limitedRho,
                        sample.rhoU,
                        sample.E,
                        epsilon
                    )
                );
        }

        const scalar meanRhoThermoEnergy =
            thermoEnergyDensity(meanRho, meanRhoU, meanE, epsilon);

        const scalar theta2 =
            limiterTheta
            (
                meanRhoThermoEnergy,
                minRhoThermoEnergy,
                min(min(epsilon, meanRho), meanRhoThermoEnergy),
                epsilon
            );

        if (theta1 >= 1.0 && theta2 >= 1.0)
        {
            continue;
        }

        scaleHighOrderModes(rhoField, cellI, theta1*theta2);
        scaleHighOrderModes(rhoUField, cellI, theta2);
        scaleHighOrderModes(EField, cellI, theta2);

        ++nLimitedCells;
    }

    reduce(nLimitedCells, sumOp<label>());

    Info<< "Lagrange-node positivity limiter applied on "
        << nLimitedCells << " cell(s)." << nl;
}
