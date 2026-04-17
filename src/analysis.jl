"""
Analysis functions specific to MatrixProjectionModels.
Matrix-level and solution-level dispatches come from StructuredPopulationCore.
"""

# --- MatrixProjectionModel forwarding dispatches ---
# Extend ProjectionModels functions with MatrixProjectionModel methods

StructuredPopulationCore.lambda(mpm::MatrixProjectionModel) = lambda(mpm.A)
StructuredPopulationCore.stable_distribution(mpm::MatrixProjectionModel) = stable_distribution(mpm.A)
StructuredPopulationCore.reproductive_value(mpm::MatrixProjectionModel) = reproductive_value(mpm.A)
StructuredPopulationCore.sensitivity(mpm::MatrixProjectionModel) = sensitivity(mpm.A)
StructuredPopulationCore.elasticity(mpm::MatrixProjectionModel) = elasticity(mpm.A)
StructuredPopulationCore.damping_ratio(mpm::MatrixProjectionModel) = damping_ratio(mpm.A)

# --- Backward-compatible alias ---

"""
    mean_matrix(sol::AbstractProjectionSolution)

Alias for `mean_kernel` — returns the mean matrix from a stochastic simulation.
"""
mean_matrix(sol::AbstractProjectionSolution) = mean_kernel(sol)
