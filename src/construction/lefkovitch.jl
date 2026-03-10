"""
Lefkovitch (stage-structured) random matrix construction.
Based on mpmsim::rand_lefko_mpm with 4 archetypes.
"""

"""
    rand_lefko_mpm(n_stages::Int, fecundity::AbstractVector;
                   archetype::Int=1, rng::AbstractRNG=Random.default_rng())

Generate a random Lefkovitch (stage-structured) matrix population model.

# Archetypes
1. All survival in stasis (diagonal): stages survive and stay
2. All survival in progression (sub-diagonal): stages survive and advance
3. Survival split between stasis and progression
4. Survival split with possible retrogression

# Arguments
- `n_stages`: Number of stages
- `fecundity`: Vector of fecundity values (length n_stages)
- `archetype`: MPM archetype (1-4)
- `rng`: Random number generator

# Returns
`MatrixProjectionModel`
"""
function rand_lefko_mpm(n_stages::Int, fecundity::AbstractVector{<:Real};
                        archetype::Int=1,
                        rng::AbstractRNG=Random.default_rng())
    1 <= archetype <= 4 || throw(ArgumentError("archetype must be 1-4"))
    length(fecundity) == n_stages || throw(ArgumentError("fecundity must have length n_stages"))
    n_stages >= 2 || throw(ArgumentError("n_stages must be >= 2"))

    T = Float64
    U = zeros(T, n_stages, n_stages)

    # Generate random survival probabilities from Uniform(0, 1)
    stage_survival = rand(rng, n_stages)

    if archetype == 1
        # All survival in stasis (diagonal)
        for i in 1:n_stages
            U[i, i] = stage_survival[i]
        end

    elseif archetype == 2
        # All survival in progression (sub-diagonal)
        for i in 1:(n_stages - 1)
            U[i+1, i] = stage_survival[i]
        end
        # Last stage: stasis
        U[n_stages, n_stages] = stage_survival[n_stages]

    elseif archetype == 3
        # Survival split between stasis and progression using Dirichlet(1,1)
        for i in 1:(n_stages - 1)
            # Split survival between staying and advancing
            split = rand(rng)
            U[i, i] = stage_survival[i] * split            # stasis
            U[i+1, i] = stage_survival[i] * (1 - split)    # progression
        end
        # Last stage: all stasis
        U[n_stages, n_stages] = stage_survival[n_stages]

    elseif archetype == 4
        # Like archetype 3 but with possible retrogression
        for i in 1:(n_stages - 1)
            if i == 1
                # First stage: split between stasis and progression
                split = rand(rng)
                U[1, 1] = stage_survival[1] * split
                U[2, 1] = stage_survival[1] * (1 - split)
            else
                # Dirichlet(1,1,1) for retrogression, stasis, progression
                d = rand(rng, Dirichlet(3, 1.0))
                U[i-1, i] = stage_survival[i] * d[1]  # retrogression
                U[i, i] = stage_survival[i] * d[2]     # stasis
                U[i+1, i] = stage_survival[i] * d[3]   # progression
            end
        end
        # Last stage: split between retrogression and stasis
        split = rand(rng)
        U[n_stages-1, n_stages] = stage_survival[n_stages] * split
        U[n_stages, n_stages] = stage_survival[n_stages] * (1 - split)
    end

    F = zeros(T, n_stages, n_stages)
    F[1, :] .= fecundity
    C = zeros(T, n_stages, n_stages)
    A = U .+ F

    stage_names = [Symbol("stage_$i") for i in 1:n_stages]
    return MatrixProjectionModel(A, U, F, C, StageClass[], stage_names)
end

"""
    rand_lefko_mpm(n_stages::Int, fecundity::Real; kwargs...)

Convenience: uniform fecundity across all stages.
"""
function rand_lefko_mpm(n_stages::Int, fecundity::Real; kwargs...)
    rand_lefko_mpm(n_stages, fill(Float64(fecundity), n_stages); kwargs...)
end
