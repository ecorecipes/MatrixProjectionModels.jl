"""
Sampling error simulation for MPMs.
Based on mpmsim::add_mpm_error.
"""

"""
    add_mpm_error(mpm::MatrixProjectionModel, sample_size::Int;
                  rng::AbstractRNG=Random.default_rng())

Add sampling error to an MPM by treating each column as a multinomial draw.
Survival transitions (U) are sampled from Binomial distributions.
Fecundity (F) is sampled from Poisson distributions.

# Arguments
- `mpm`: Input MatrixProjectionModel
- `sample_size`: Number of individuals sampled per stage
- `rng`: Random number generator

# Returns
New `MatrixProjectionModel` with sampling error applied.
"""
function add_mpm_error(mpm::MatrixProjectionModel, sample_size::Int;
                       rng::AbstractRNG=Random.default_rng())
    sample_size > 0 || throw(ArgumentError("sample_size must be positive"))
    n = n_stages(mpm)

    U_new = zeros(Float64, n, n)
    F_new = zeros(Float64, n, n)
    C_new = zeros(Float64, n, n)

    # Survival transitions: each element U[i,j] sampled via Binomial
    # Each entry is an independent Binomial trial: out of sample_size individuals,
    # how many transition from j to i?
    for j in 1:n
        for i in 1:n
            if mpm.U[i, j] > 0
                p = min(mpm.U[i, j], 1.0)
                count = rand(rng, Binomial(sample_size, p))
                U_new[i, j] = count / sample_size
            end
        end
    end

    # Fecundity: Poisson sampling
    for j in 1:n
        for i in 1:n
            if mpm.F[i, j] > 0
                total_offspring = rand(rng, Poisson(mpm.F[i, j] * sample_size))
                F_new[i, j] = total_offspring / sample_size
            end
        end
    end

    # Clonal: same as fecundity
    for j in 1:n
        for i in 1:n
            if mpm.C[i, j] > 0
                total_offspring = rand(rng, Poisson(mpm.C[i, j] * sample_size))
                C_new[i, j] = total_offspring / sample_size
            end
        end
    end

    return MatrixProjectionModel(U_new, F_new, C_new;
                                 stages=mpm.stages, stage_names=mpm.stage_names)
end

"""
    add_mpm_error(A::AbstractMatrix, sample_size::Int;
                  rng::AbstractRNG=Random.default_rng())

Add sampling error to a raw projection matrix (treats all as survival transitions).
"""
function add_mpm_error(A::AbstractMatrix, sample_size::Int;
                       rng::AbstractRNG=Random.default_rng())
    mpm = MatrixProjectionModel(A)
    return add_mpm_error(mpm, sample_size; rng=rng)
end
