// Create reusable DG PDE terms once at solver setup.
//
// Each explicit stage supplies the active cell fields through add(U, F), then
// updates the numerical face flux with updateFlux(cellI).
const bool selfDiffusion = thermo->selfDiffusion();
const bool artificialDiffusion = artificialDiffusionNonSmooth.enabled();

dgConvectiveTerm<scalar, vector> massConvTerm
(
    "massConvTerm",
    dgFluxSolver::equationType::massTransport,
    dgFluxSolver::fluxType::convective,
    dgMesh,
    fluxSolverManager
);

std::optional<dgDiffusiveTerm<scalar, vector>> selfDiffusionMassTerm;

if (selfDiffusion)
{
    selfDiffusionMassTerm.emplace
    (
        "selfDiffusionMassTerm",
        dgFluxSolver::equationType::massTransport,
        dgFluxSolver::fluxType::diffusive,
        dgMesh,
        fluxSolverManager
    );
}

std::optional<dgDiffusiveTerm<scalar, vector>> artificialDiffusionMassTerm;

if (artificialDiffusion)
{
    artificialDiffusionMassTerm.emplace
    (
        "artificialDiffusionMassTerm",
        dgFluxSolver::equationType::massTransport,
        dgFluxSolver::fluxType::diffusive,
        dgMesh,
        fluxSolverManager
    );
}

dgConvectiveTerm<vector, tensor> momentumConvTerm
(
    "momentumConvTerm",
    dgFluxSolver::equationType::momentumTransport,
    dgFluxSolver::fluxType::convective,
    dgMesh,
    fluxSolverManager
);

dgDiffusiveTerm<vector, tensor> momentumDiffTerm
(
    "momentumDiffTerm",
    dgFluxSolver::equationType::momentumTransport,
    dgFluxSolver::fluxType::diffusive,
    dgMesh,
    fluxSolverManager
);

std::optional<dgDiffusiveTerm<vector, tensor>>
    artificialDiffusionMomentumTerm;

if (artificialDiffusion)
{
    artificialDiffusionMomentumTerm.emplace
    (
        "artificialDiffusionMomentumTerm",
        dgFluxSolver::equationType::momentumTransport,
        dgFluxSolver::fluxType::diffusive,
        dgMesh,
        fluxSolverManager
    );
}

dgConvectiveTerm<scalar, vector> energyConvTerm
(
    "energyConvTerm",
    dgFluxSolver::equationType::energyTransport,
    dgFluxSolver::fluxType::convective,
    dgMesh,
    fluxSolverManager
);

dgDiffusiveTerm<scalar, vector> energyDiffTerm
(
    "energyDiffTerm",
    dgFluxSolver::equationType::energyTransport,
    dgFluxSolver::fluxType::diffusive,
    dgMesh,
    fluxSolverManager
);

dgSourceTerm<scalar> kineticMassDiffCorrectionTerm
(
    "kineticMassDiffCorrectionTerm",
    dgMesh
);

std::optional<dgDiffusiveTerm<scalar, vector>>
    artificialDiffusionEnergyTerm;

if (artificialDiffusion)
{
    artificialDiffusionEnergyTerm.emplace
    (
        "artificialDiffusionEnergyTerm",
        dgFluxSolver::equationType::energyTransport,
        dgFluxSolver::fluxType::diffusive,
        dgMesh,
        fluxSolverManager
    );
}
