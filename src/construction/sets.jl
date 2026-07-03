"""
Batch generation of MPM sets.
Based on mpmsim::rand_leslie_set and mpmsim::rand_lefko_set.
"""

"""
    rand_leslie_set(n::Int; mortality_model, fecundity_model,
                    mortality_params_dist=nothing, fecundity_params_dist=nothing,
                    n_stages=nothing, truncate=0.01, rng=Random.default_rng())

Generate a set of n random Leslie MPMs by sampling parameters from distributions.

If no parameter distributions are provided, generates n identical MPMs from
the given mortality and fecundity models.
"""
_sample_param(rng::AbstractRNG, value) = value
_sample_param(rng::AbstractRNG, dist::Distribution) = rand(rng, dist)

function _resample_model(model, params_dist, rng::AbstractRNG)
    params_dist === nothing && return model
    names = fieldnames(typeof(model))
    values = ntuple(length(names)) do i
        name = names[i]
        current = getfield(model, name)
        source = if params_dist isa NamedTuple
            haskey(params_dist, name) ? getproperty(params_dist, name) : current
        elseif params_dist isa Tuple || params_dist isa AbstractVector
            length(params_dist) == length(names) ||
                throw(ArgumentError("parameter collections must have length $(length(names)) for $(typeof(model))"))
            params_dist[i]
        else
            length(names) == 1 ||
                throw(ArgumentError("scalar/distribution parameter specifications are only supported for single-parameter models"))
            params_dist
        end
        _sample_param(rng, source)
    end
    return typeof(model)(values...)
end

function rand_leslie_set(n::Int;
                         mortality_model::AbstractMortalityModel,
                         fecundity_model::AbstractFecundityModel,
                         mortality_params_dist=nothing,
                         fecundity_params_dist=nothing,
                         n_stages::Union{Nothing,Int}=nothing,
                         truncate::Real=0.01,
                         rng::AbstractRNG=Random.default_rng())
    mpms = Vector{MatrixProjectionModel{Float64}}(undef, n)
    for i in 1:n
        sampled_mortality = _resample_model(mortality_model, mortality_params_dist, rng)
        sampled_fecundity = _resample_model(fecundity_model, fecundity_params_dist, rng)
        mpms[i] = make_leslie_mpm(sampled_mortality, sampled_fecundity;
                                  n_stages=n_stages, truncate=truncate)
    end
    return mpms
end

"""
    rand_lefko_set(n::Int; n_stages::Int, fecundity, archetype::Int=1,
                   constraint=nothing, rng=Random.default_rng())

Generate a set of n random Lefkovitch MPMs.

# Arguments
- `n`: Number of MPMs to generate
- `n_stages`: Number of stages
- `fecundity`: Scalar or vector of fecundity values
- `archetype`: MPM archetype (1-4)
- `constraint`: Optional function `f(mpm) -> Bool` to filter; if provided,
  keeps generating until n valid MPMs are found.
- `rng`: Random number generator
"""
function rand_lefko_set(n::Int;
                        n_stages::Int,
                        fecundity,
                        archetype::Int=1,
                        constraint=nothing,
                        rng::AbstractRNG=Random.default_rng())
    fec_vec = fecundity isa Real ? fill(Float64(fecundity), n_stages) : fecundity

    mpms = Vector{MatrixProjectionModel{Float64}}()
    max_attempts = n * 100  # Prevent infinite loops

    attempts = 0
    while length(mpms) < n && attempts < max_attempts
        attempts += 1
        mpm = rand_lefko_mpm(n_stages, fec_vec; archetype=archetype, rng=rng)
        if constraint === nothing || constraint(mpm)
            push!(mpms, mpm)
        end
    end

    length(mpms) < n &&
        @warn "Only generated $(length(mpms)) of $n requested MPMs after $max_attempts attempts"

    return mpms
end
