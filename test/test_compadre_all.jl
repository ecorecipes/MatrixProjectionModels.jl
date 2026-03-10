"""
Test MatrixProjectionModels.jl against all COMPADRE and COMADRE matrices.

Prerequisites:
  1. Run the R export scripts to create compadre_plain.rds and comadre_plain.rds in a temp directory
  2. Or place pre-existing JLD2 caches in a temp directory

Run: julia --project test/test_compadre_all.jl
"""

using CodecBzip2, RData, JLD2
using MatrixProjectionModels
using LinearAlgebra, Statistics, Printf

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

    isfile(rds_path) || error("No data for $db_name. Export from R first.")
    @info "Loading $db_name from RDS and caching to JLD2"
    d = RData.load(rds_path)
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

function test_database(db_name)
    println("\n" * "="^70)
    println("Testing $(uppercase(db_name))")
    println("="^70)

    db = load_db(db_name)
    n = length(db.species)
    println("Records: $n")

    n_tested = 0; n_lambda_exact = 0; n_lambda_close = 0; n_lambda_ok = 0
    n_build_ok = 0; n_build_fail = 0; n_analysis_ok = 0; n_analysis_fail = 0
    lambda_errors = Float64[]
    t0 = time()

    for i in 1:n
        db.dims[i] == 0 && continue
        r_lam = db.lambda_R[i]
        isnan(r_lam) && continue
        n_tested += 1

        try
            mpm = MatrixProjectionModel(db.A[i], db.U[i], db.F[i], db.C[i])
            n_build_ok += 1

            jl_lam = lambda(mpm)
            if isfinite(r_lam) && isfinite(jl_lam)
                rel_err = abs(jl_lam - r_lam) / max(abs(r_lam), 1e-10)
                push!(lambda_errors, rel_err)
                rel_err < 1e-10 && (n_lambda_exact += 1)
                rel_err < 0.01 && (n_lambda_close += 1)
                rel_err < 0.05 && (n_lambda_ok += 1)
            end

            if n_tested <= 300 || n_tested % 20 == 0
                try
                    sensitivity(mpm); elasticity(mpm); stable_distribution(mpm)
                    n_analysis_ok += 1
                catch
                    n_analysis_fail += 1
                end
            end
        catch
            n_build_fail += 1
        end

        n_tested % 2000 == 0 && @printf("  %d/%d (%.1fs)...\n", n_tested, n, time()-t0)
    end

    elapsed = time() - t0
    println("\n--- $(uppercase(db_name)) Results ---")
    @printf("Models tested:      %6d / %d\n", n_tested, n)
    @printf("Build success:      %6d / %d  (%.1f%%)\n", n_build_ok, n_tested, 100n_build_ok/max(n_tested,1))
    @printf("Lambda exact:       %6d / %d  (%.1f%%)\n", n_lambda_exact, n_tested, 100n_lambda_exact/max(n_tested,1))
    @printf("Lambda <1%% error:   %6d / %d  (%.1f%%)\n", n_lambda_close, n_tested, 100n_lambda_close/max(n_tested,1))
    @printf("Lambda <5%% error:   %6d / %d  (%.1f%%)\n", n_lambda_ok, n_tested, 100n_lambda_ok/max(n_tested,1))
    if !isempty(lambda_errors)
        @printf("Relative error:     mean=%.2e  median=%.2e  max=%.2e\n",
                mean(lambda_errors), median(lambda_errors), maximum(lambda_errors))
    end
    @printf("Analysis tests:     %d pass, %d fail\n", n_analysis_ok, n_analysis_fail)
    @printf("Time: %.1f seconds\n", elapsed)

    return (n_tested=n_tested, n_build_ok=n_build_ok,
            n_lambda_close=n_lambda_close, n_lambda_ok=n_lambda_ok)
end

function test_life_history(db_name; n_sample=500)
    println("\n--- Life History: $(uppercase(db_name)) (n=$n_sample) ---")
    db = load_db(db_name)
    valid = [i for i in 1:length(db.species)
             if !isnan(db.lambda_R[i]) && 0.5 < db.lambda_R[i] < 2.0 && 2 <= db.dims[i] <= 15]
    ids = valid[1:min(n_sample, length(valid))]

    n_ok = 0; n_fail = 0
    traits = Dict{String, Vector{Float64}}()

    for i in ids
        try
            U = db.U[i]; Fm = db.F[i]
            le = life_expect_mean(U)
            R0 = net_repro_rate(U, Fm; method=:fundamental)
            Tg = gen_time(U, Fm; method=:cohort)
            lo = longevity(U)
            surv = vr_survival(U)

            for (k, v) in [("life_expect", le), ("R0", R0), ("gen_time", Tg),
                           ("longevity", Float64(lo)), ("survival", surv)]
                isfinite(v) && push!(get!(traits, k, Float64[]), v)
            end
            n_ok += 1
        catch
            n_fail += 1
        end
    end

    @printf("Success: %d / %d (%.1f%%)\n", n_ok, length(ids), 100n_ok/max(length(ids),1))
    for (k, v) in sort(collect(traits))
        @printf("  %-15s  n=%5d  mean=%8.3f  median=%8.3f  [%.3f, %.3f]\n",
                k, length(v), mean(v), median(v), minimum(v), maximum(v))
    end
end

# === Main ===
println("MatrixProjectionModels.jl — Full COMPADRE/COMADRE Validation")
println("="^70)

results = Dict{String, Any}()
for db in ["compadre", "comadre"]
    try
        results[db] = test_database(db)
        test_life_history(db)
    catch e
        println("ERROR ($db): $e")
    end
end

println("\n" * "="^70)
println("SUMMARY")
println("="^70)
for (db, r) in sort(collect(results))
    @printf("%-10s  %d tested, %d build OK, %d/%d (%.1f%%) lambda <1%%\n",
            uppercase(db), r.n_tested, r.n_build_ok, r.n_lambda_close, r.n_tested,
            100r.n_lambda_close/max(r.n_tested,1))
end
