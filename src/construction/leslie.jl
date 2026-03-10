"""
Leslie matrix construction from survival and fecundity schedules.
Based on mpmsim::make_leslie_mpm.
"""

"""
    make_leslie_mpm(survival::AbstractVector, fecundity::AbstractVector; split=true)

Construct a Leslie matrix population model from survival probabilities and
fecundity values.

# Arguments
- `survival`: Vector of survival probabilities for each age class (length n-1 or n)
- `fecundity`: Vector of fecundity values for each age class (length n)
- `split`: If true, decompose into U (survival) and F (fecundity). Default: true.

# Returns
`MatrixProjectionModel` with Leslie structure.
"""
function make_leslie_mpm(survival::AbstractVector{<:Real},
                         fecundity::AbstractVector{<:Real};
                         split::Bool=true)
    n = length(fecundity)
    ns = length(survival)
    (ns == n - 1 || ns == n) ||
        throw(ArgumentError("survival must have length $(n-1) or $n, got $ns"))

    T = promote_type(eltype(survival), eltype(fecundity))
    A = zeros(T, n, n)
    U = zeros(T, n, n)
    F = zeros(T, n, n)

    # Fecundity in first row
    F[1, :] .= fecundity

    # Survival on sub-diagonal
    for i in 1:min(ns, n-1)
        U[i+1, i] = survival[i]
    end

    A .= U .+ F
    C = zeros(T, n, n)
    stage_names = [Symbol("age_$i") for i in 1:n]

    return MatrixProjectionModel(A, U, F, C, StageClass[], stage_names)
end

"""
    make_leslie_mpm(mortality_model::AbstractMortalityModel,
                    fecundity_model::AbstractFecundityModel;
                    n_stages=nothing, truncate=0.01)

Construct a Leslie MPM from mortality and fecundity models.
Number of stages is determined by truncation of survivorship (lx < truncate).
"""
function make_leslie_mpm(mortality_model::AbstractMortalityModel,
                         fecundity_model::AbstractFecundityModel;
                         n_stages::Union{Nothing,Int}=nothing,
                         truncate::Real=0.01)
    surv = model_survival(mortality_model; truncate=truncate)
    fec = model_fecundity(fecundity_model; ages=surv.x)

    n = n_stages !== nothing ? n_stages : length(surv.x)
    n = min(n, length(surv.x))

    survival = surv.px[1:n-1]
    fecundity = fec.fx[1:n]

    return make_leslie_mpm(survival, fecundity)
end
