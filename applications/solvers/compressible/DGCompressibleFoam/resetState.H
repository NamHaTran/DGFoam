// Some flux solvers cache face-local intermediate states to avoid recomputing
// them multiple times within one assembly pass. Reset that cache whenever the
// stage input state changes. This keeps each stage assembly self-consistent.
fluxSolverManager.reset();
