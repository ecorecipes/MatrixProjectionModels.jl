"""
Maturity/reproductive onset calculations.
Based on Rage::mature_prob, mature_age, mature_distrib.
"""

"""
    mature_prob(matU::AbstractMatrix, matR::AbstractMatrix; start::Int=1)

Probability of reaching reproductive maturity starting from a given stage.
matR is the reproductive matrix (F or F+C).
"""
function mature_prob(matU::AbstractMatrix, matR::AbstractMatrix; start::Int=1)
    n = size(matU, 1)
    _validate_start(start, n)

    # Identify reproductive stages
    repro = repro_stages(matR)
    any(repro) || return 0.0

    # Probability of reaching any reproductive stage
    N = _fundamental_matrix(matU)
    N === nothing && return 0.0

    # Sum of expected time in reproductive stages, weighted by entry probability
    # More precisely: prob of ever entering a repro stage
    # Use absorbing Markov chain: make repro stages absorbing
    matU_abs = copy(Float64.(matU))
    for j in 1:n
        if repro[j]
            matU_abs[:, j] .= 0.0
            matU_abs[j, j] = 1.0
        end
    end

    # Iterate until convergence
    state = zeros(Float64, n)
    state[start] = 1.0

    for _ in 1:10000
        new_state = matU_abs * state
        total_in_repro = sum(new_state[j] for j in 1:n if repro[j])
        if total_in_repro > 1 - 1e-10 || sum(new_state) < 1e-15
            return min(total_in_repro, 1.0)
        end
        state = new_state
    end

    return sum(state[j] for j in 1:n if repro[j])
end

"""
    mature_age(matU::AbstractMatrix, matR::AbstractMatrix; start::Int=1)

Mean age at reproductive maturity.
"""
function mature_age(matU::AbstractMatrix, matR::AbstractMatrix; start::Int=1)
    n = size(matU, 1)
    _validate_start(start, n)

    repro = repro_stages(matR)
    any(repro) || return Inf

    # Non-reproductive stages
    non_repro = .!repro

    # Extract sub-matrix for non-reproductive stages
    nr_idx = findall(non_repro)
    if isempty(nr_idx)
        return 0.0  # Already reproductive
    end

    if non_repro[start]
        # Fundamental matrix of non-reproductive sub-matrix
        U_nr = matU[nr_idx, nr_idx]
        N_nr = _fundamental_matrix(U_nr)
        N_nr === nothing && return Inf

        # Find position of start in nr_idx
        start_pos = findfirst(==(start), nr_idx)
        start_pos === nothing && return 0.0

        return sum(N_nr[start_pos, :])
    else
        return 0.0  # Start stage is already reproductive
    end
end

"""
    mature_distrib(matU::AbstractMatrix; repro_stages::AbstractVector{Bool})

Distribution of individuals across stages at first reproduction.
"""
function mature_distrib(matU::AbstractMatrix; repro_stages::AbstractVector{Bool})
    n = size(matU, 1)

    non_repro = .!repro_stages
    nr_idx = findall(non_repro)
    r_idx = findall(repro_stages)

    isempty(nr_idx) && return ones(n) / n
    isempty(r_idx) && return zeros(n)

    U_nr = matU[nr_idx, nr_idx]
    N_nr = _fundamental_matrix(U_nr)
    N_nr === nothing && return zeros(n)

    # Transition from non-repro to repro
    T_rn = matU[r_idx, nr_idx]

    # Distribution of first entry into reproductive stages
    # d = T_rn * N_nr * initial
    dist = zeros(n)
    # Assume starting from first non-repro stage
    initial = zeros(length(nr_idx))
    initial[1] = 1.0

    entry = T_rn * N_nr * initial
    for (k, j) in enumerate(r_idx)
        dist[j] = max(entry[k], 0.0)
    end

    s = sum(dist)
    s > 0 && (dist ./= s)
    return dist
end
