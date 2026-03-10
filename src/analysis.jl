"""
Analysis functions specific to MatrixProjectionModels.
Matrix-level and solution-level dispatches come from ProjectionModels.
"""

# --- MatrixProjectionModel forwarding dispatches ---
# Extend ProjectionModels functions with MatrixProjectionModel methods

ProjectionModels.lambda(mpm::MatrixProjectionModel) = lambda(mpm.A)
ProjectionModels.stable_distribution(mpm::MatrixProjectionModel) = stable_distribution(mpm.A)
ProjectionModels.reproductive_value(mpm::MatrixProjectionModel) = reproductive_value(mpm.A)
ProjectionModels.sensitivity(mpm::MatrixProjectionModel) = sensitivity(mpm.A)
ProjectionModels.elasticity(mpm::MatrixProjectionModel) = elasticity(mpm.A)
ProjectionModels.damping_ratio(mpm::MatrixProjectionModel) = damping_ratio(mpm.A)

# --- Backward-compatible alias ---

"""
    mean_matrix(sol::AbstractProjectionSolution)

Alias for `mean_kernel` — returns the mean matrix from a stochastic simulation.
"""
mean_matrix(sol::AbstractProjectionSolution) = mean_kernel(sol)
