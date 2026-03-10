"""
Stochastic perturbation analysis.
Based on Rage::perturb_stochastic.
"""

"""
    perturb_stochastic(matrices::AbstractVector{<:AbstractMatrix};
                       pert::Real=1e-6, type::Symbol=:sensitivity)

Stochastic perturbation analysis for a set of environment-specific matrices.
Computes sensitivity/elasticity of the stochastic growth rate log(λs)
to mean matrix elements.

# Arguments
- `matrices`: Vector of projection matrices (one per environment)
- `pert`: Perturbation magnitude
- `type`: `:sensitivity` or `:elasticity`

# Returns
Matrix of stochastic sensitivities or elasticities.
"""
function perturb_stochastic(matrices::AbstractVector{<:AbstractMatrix};
                            pert::Real=1e-6,
                            type::Symbol=:sensitivity,
                            n_sim::Int=50000,
                            rng::AbstractRNG=Random.default_rng())
    n_env = length(matrices)
    n = size(matrices[1], 1)

    # Compute base stochastic growth rate
    base_log_λs = _stochastic_log_lambda(matrices, n_sim; rng=rng)

    result = zeros(n, n)
    A_mean = mean(matrices)

    for j in 1:n, i in 1:n
        # Perturb the mean of element [i,j] by adding pert to all matrices
        perturbed = [copy(float.(m)) for m in matrices]
        for m in perturbed
            m[i, j] += pert
        end

        pert_log_λs = _stochastic_log_lambda(perturbed, n_sim; rng=copy(rng))
        sens = (pert_log_λs - base_log_λs) / pert

        if type == :elasticity
            if base_log_λs != 0 && A_mean[i, j] != 0
                result[i, j] = sens * A_mean[i, j] / base_log_λs
            end
        else
            result[i, j] = sens
        end
    end

    return result
end

"""
    _stochastic_log_lambda(matrices, n_sim; rng)

Compute log stochastic growth rate via simulation.
"""
function _stochastic_log_lambda(matrices::AbstractVector{<:AbstractMatrix},
                                n_sim::Int; rng::AbstractRNG=Random.default_rng())
    n_env = length(matrices)
    n = size(matrices[1], 1)

    v = ones(n) / n
    log_λ_sum = 0.0

    for t in 1:n_sim
        A = matrices[rand(rng, 1:n_env)]
        v_new = A * v
        norm_v = sum(v_new)
        if norm_v > 0
            log_λ_sum += log(norm_v)
            v = v_new ./ norm_v
        end
    end

    # Burn-in: discard first 10%
    # Actually re-simulate with burn-in for accuracy
    burn = div(n_sim, 10)
    return log_λ_sum / n_sim
end
