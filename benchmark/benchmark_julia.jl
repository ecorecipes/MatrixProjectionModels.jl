"""
Benchmark Julia MatrixProjectionModels.jl analysis.
Uses BenchmarkTools.jl for precise timing.
Data: COMPADRE/COMADRE pre-exported as plain RDS, cached as JLD2.

Run: julia --project benchmark/benchmark_julia.jl
"""

using BenchmarkTools
using JLD2
using MatrixProjectionModels
using LinearAlgebra, Statistics, Printf

# RData is a weak dep; only needed if JLD2 cache doesn't exist
const HAS_RDATA = try
    @eval using RData
    true
catch
    false
end

const CACHE_DIR = joinpath(tempdir(), "compadre_data")

function to_matrix(raw, dim)
    raw isa Number && return fill(Float64(raw), 1, 1)
    isempty(raw) && return zeros(Float64, dim, dim)
    v = Float64.(raw)
    len = length(v)
    if len == dim * dim
        return reshape(v, dim, dim)
    end
    actual_dim = isqrt(len)
    if actual_dim * actual_dim == len
        return reshape(v, actual_dim, actual_dim)
    end
    return zeros(Float64, dim, dim)
end

function load_db(db_name)
    jld2_path = joinpath(CACHE_DIR, "$(db_name).jld2")
    rds_path = joinpath(tempdir(), "$(db_name)_plain.rds")

    if isfile(jld2_path)
        @info "Loading $db_name from JLD2 cache"
        return jldopen(jld2_path, "r") do f
            (species=f["species"], dims=f["dims"], lambda_R=f["lambda_R"],
             A=f["A"], U=f["U"], F=f["F"], C=f["C"])
        end
    end

    isfile(rds_path) || error("No data at $rds_path or $jld2_path")
    HAS_RDATA || error("RData.jl not available. Run the test_compadre_all.jl script first to create JLD2 caches.")
    @info "Loading $db_name from RDS and caching to JLD2"
    d = Main.RData.load(rds_path)
    species = String.(d["species"])
    dims = Int.(d["dim"])
    lambda_R = Float64[ismissing(x) || (x isa Number && isnan(x)) ? NaN : Float64(x) for x in d["lambda"]]
    n = length(species)

    A = Matrix{Float64}[]
    U = Matrix{Float64}[]
    Fm = Matrix{Float64}[]
    Cm = Matrix{Float64}[]
    for i in 1:n
        dim = dims[i] == 0 ? 1 : dims[i]
        raw_a = d["matA"][i]
        if !(raw_a isa Number) && !isempty(raw_a)
            actual_dim = isqrt(length(raw_a))
            if actual_dim * actual_dim == length(raw_a)
                dim = actual_dim
            end
        end
        dims[i] = dim
        push!(A, to_matrix(d["matA"][i], dim))
        push!(U, to_matrix(d["matU"][i], dim))
        push!(Fm, to_matrix(d["matF"][i], dim))
        push!(Cm, to_matrix(d["matC"][i], dim))
    end

    mkpath(CACHE_DIR)
    jldsave(jld2_path; species, dims, lambda_R, A, U, F=Fm, C=Cm)
    return (species=species, dims=dims, lambda_R=lambda_R, A=A, U=U, F=Fm, C=Cm)
end

function select_models(db; n_sample=500)
    valid = [i for i in 1:length(db.species)
             if !isnan(db.lambda_R[i]) && 0.1 < db.lambda_R[i] < 10.0 &&
                2 <= db.dims[i] <= 20]
    return valid[1:min(n_sample, length(valid))]
end

# --- Full analysis pipeline (matching R benchmark) ---
function full_analysis(A, U, Fm)
    lam = lambda(A)
    sens = sensitivity(A)
    elas = elasticity(A)
    le = life_expect_mean(U)
    R0 = net_repro_rate(U, Fm; method=:fundamental)
    Tg = gen_time(U, Fm; method=:R0)
    lo = longevity(U)
    vrs = vr_survival(U)
    (lambda=lam, le=le, R0=R0, gen_time=Tg, longevity=lo, vr_surv=vrs)
end

function run_benchmark(db_name)
    println("\n" * "="^60)
    println("Benchmarking $(uppercase(db_name))")
    println("="^60)

    db = load_db(db_name)
    ids = select_models(db)
    println("Selected $(length(ids)) models for benchmarking")

    # Pre-extract matrices for the selected models
    As = [db.A[i] for i in ids]
    Us = [db.U[i] for i in ids]
    Fs = [db.F[i] for i in ids]

    # Warmup
    println("Warming up...")
    for i in 1:min(10, length(ids))
        try; full_analysis(As[i], Us[i], Fs[i]); catch; end
    end

    # --- Batch benchmarks (all models) ---
    println("\n--- Batch benchmarks (all $(length(ids)) models, BenchmarkTools) ---")

    # Lambda
    b = @benchmark begin
        for i in 1:length($As)
            try; lambda($As[i]); catch; end
        end
    end samples=10 evals=1
    t_lambda = median(b.times) / 1e9
    @printf("  lambda:         %.4f s  (%d models, %.1f us/model)\n",
            t_lambda, length(ids), t_lambda / length(ids) * 1e6)

    # Sensitivity
    b = @benchmark begin
        for i in 1:length($As)
            try; sensitivity($As[i]); catch; end
        end
    end samples=10 evals=1
    t_sens = median(b.times) / 1e9
    @printf("  sensitivity:    %.4f s  (%d models, %.1f us/model)\n",
            t_sens, length(ids), t_sens / length(ids) * 1e6)

    # Elasticity
    b = @benchmark begin
        for i in 1:length($As)
            try; elasticity($As[i]); catch; end
        end
    end samples=10 evals=1
    t_elas = median(b.times) / 1e9
    @printf("  elasticity:     %.4f s  (%d models, %.1f us/model)\n",
            t_elas, length(ids), t_elas / length(ids) * 1e6)

    # Life expectancy
    b = @benchmark begin
        for i in 1:length($Us)
            try; life_expect_mean($Us[i]); catch; end
        end
    end samples=10 evals=1
    t_le = median(b.times) / 1e9
    @printf("  life_expect:    %.4f s  (%d models, %.1f us/model)\n",
            t_le, length(ids), t_le / length(ids) * 1e6)

    # Net reproductive rate
    b = @benchmark begin
        for i in 1:length($Us)
            try; net_repro_rate($Us[i], $Fs[i]; method=:fundamental); catch; end
        end
    end samples=10 evals=1
    t_R0 = median(b.times) / 1e9
    @printf("  net_repro_rate: %.4f s  (%d models, %.1f us/model)\n",
            t_R0, length(ids), t_R0 / length(ids) * 1e6)

    # Generation time
    b = @benchmark begin
        for i in 1:length($Us)
            try; gen_time($Us[i], $Fs[i]; method=:R0); catch; end
        end
    end samples=10 evals=1
    t_gen = median(b.times) / 1e9
    @printf("  gen_time:       %.4f s  (%d models, %.1f us/model)\n",
            t_gen, length(ids), t_gen / length(ids) * 1e6)

    # Longevity
    b = @benchmark begin
        for i in 1:length($Us)
            try; longevity($Us[i]); catch; end
        end
    end samples=10 evals=1
    t_long = median(b.times) / 1e9
    @printf("  longevity:      %.4f s  (%d models, %.1f us/model)\n",
            t_long, length(ids), t_long / length(ids) * 1e6)

    # Full pipeline
    b = @benchmark begin
        for i in 1:length($As)
            try; full_analysis($As[i], $Us[i], $Fs[i]); catch; end
        end
    end samples=10 evals=1
    t_full = median(b.times) / 1e9
    @printf("  full_pipeline:  %.4f s  (%d models, %.1f us/model)\n",
            t_full, length(ids), t_full / length(ids) * 1e6)

    # --- Single matrix benchmarks (representative sizes) ---
    println("\n--- Single matrix benchmarks (BenchmarkTools @btime) ---")
    sizes = [3, 5, 10]
    for sz in sizes
        idx = findfirst(i -> db.dims[ids[i]] == sz, 1:length(ids))
        idx === nothing && continue
        i = idx
        A = As[i]; U = Us[i]; Fm = Fs[i]

        b_lam = @benchmark lambda($A)
        @printf("  lambda %dx%d:      median=%.1f us  (min=%.1f us)\n",
                sz, sz, median(b_lam.times)/1e3, minimum(b_lam.times)/1e3)

        b_sens = @benchmark sensitivity($A)
        @printf("  sensitivity %dx%d: median=%.1f us  (min=%.1f us)\n",
                sz, sz, median(b_sens.times)/1e3, minimum(b_sens.times)/1e3)

        b_full = @benchmark full_analysis($A, $U, $Fm)
        @printf("  full %dx%d:        median=%.1f us  (min=%.1f us)\n",
                sz, sz, median(b_full.times)/1e3, minimum(b_full.times)/1e3)
    end

    return (n=length(ids), t_lambda=t_lambda, t_sens=t_sens, t_elas=t_elas,
            t_le=t_le, t_R0=t_R0, t_gen=t_gen, t_long=t_long, t_full=t_full)
end

# --- Main ---
println("MatrixProjectionModels.jl — Benchmark")
println("="^60)

results = Dict{String, Any}()
for db_name in ["compadre", "comadre"]
    try
        results[db_name] = run_benchmark(db_name)
    catch e
        println("ERROR ($db_name): $e")
        println(sprint(showerror, e, catch_backtrace()))
    end
end

println("\n" * "="^60)
println("SUMMARY")
println("="^60)
for (db_name, r) in sort(collect(results))
    @printf("%-10s  %d models  full_pipeline=%.4fs (%.1f us/model)\n",
            uppercase(db_name), r.n, r.t_full, r.t_full/r.n * 1e6)
end

# Save for comparison
const RESULTS_PATH = joinpath(tempdir(), "benchmark_julia_mpm_results.jld2")
jldsave(RESULTS_PATH; results)
println("\nResults saved to $RESULTS_PATH")
