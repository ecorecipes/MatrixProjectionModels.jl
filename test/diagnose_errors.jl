"""
Diagnose models with >1% lambda error vs R.

Key hypothesis: R computes lambda from stored A matrix, while Julia
reconstructs A = U + F + C. If A ≠ U + F + C (rounding, unaccounted
components), lambdas will differ.

Run: julia --project test/diagnose_errors.jl
"""

using JLD2, RData
using MatrixProjectionModels
using LinearAlgebra, Statistics, Printf

const CACHE_DIR = joinpath(tempdir(), "compadre_data")

function to_matrix(raw, dim)
    raw isa Number && return fill(Float64(raw), 1, 1)
    isempty(raw) && return zeros(Float64, dim, dim)
    v = Float64.(raw)
    len = length(v)
    # If vector length matches dim²
    if len == dim * dim
        return reshape(v, dim, dim)
    end
    # Try to infer actual dimension from vector length
    actual_dim = isqrt(len)
    if actual_dim * actual_dim == len
        return reshape(v, actual_dim, actual_dim)
    end
    # Fallback: zero matrix of claimed dim
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

    isfile(rds_path) || error("No data for $db_name at $rds_path or $jld2_path")
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
        # Determine actual dimension from matA vector length
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

function diagnose(db_name)
    println("\n" * "="^70)
    println("DIAGNOSING $(uppercase(db_name))")
    println("="^70)

    db = load_db(db_name)
    n = length(db.species)

    # Categories
    cat_A_mismatch = Int[]       # A ≠ U + F + C
    cat_A_match_err = Int[]      # A = U+F+C but still error (shouldn't happen)
    cat_negative = Int[]         # Matrices have negative entries
    cat_zero_lambda = Int[]      # R lambda ≈ 0
    cat_huge_lambda = Int[]      # R lambda very large
    cat_all_zero = Int[]         # A is all zeros
    cat_numerical = Int[]        # Small numerical differences
    errors_from_A = Float64[]    # Error when computing from stored A directly
    errors_from_UFC = Float64[]  # Error when computing from U+F+C

    n_over1pct = 0

    for i in 1:n
        db.dims[i] == 0 && continue
        r_lam = db.lambda_R[i]
        isnan(r_lam) && continue
        !isfinite(r_lam) && continue

        Ai = db.A[i]
        Ui = db.U[i]
        Fi = db.F[i]
        Ci = db.C[i]
        UFC = Ui .+ Fi .+ Ci

        # Julia lambda from U+F+C (what MatrixProjectionModel does)
        jl_lam_ufc = try
            maximum(real.(eigen(UFC).values))
        catch
            NaN
        end

        # Julia lambda from stored A directly
        jl_lam_A = try
            maximum(real.(eigen(Ai).values))
        catch
            NaN
        end

        !isfinite(jl_lam_ufc) && continue

        rel_err_ufc = abs(jl_lam_ufc - r_lam) / max(abs(r_lam), 1e-10)
        rel_err_A = isfinite(jl_lam_A) ? abs(jl_lam_A - r_lam) / max(abs(r_lam), 1e-10) : NaN

        push!(errors_from_A, isfinite(rel_err_A) ? rel_err_A : NaN)
        push!(errors_from_UFC, rel_err_ufc)

        rel_err_ufc < 0.01 && continue
        n_over1pct += 1

        # Check A vs U+F+C mismatch
        a_diff = maximum(abs.(Ai .- UFC))
        has_negative = any(x -> x < -1e-10, Ai) || any(x -> x < -1e-10, UFC)

        if a_diff > 1e-10
            push!(cat_A_mismatch, i)
        else
            push!(cat_A_match_err, i)
        end

        if has_negative
            push!(cat_negative, i)
        end

        if abs(r_lam) < 1e-6
            push!(cat_zero_lambda, i)
        end

        if abs(r_lam) > 100
            push!(cat_huge_lambda, i)
        end

        if maximum(abs.(Ai)) < 1e-10
            push!(cat_all_zero, i)
        end

        # Check if using stored A fixes the error
        if isfinite(rel_err_A) && rel_err_A < 0.01
            push!(cat_numerical, i)
        end
    end

    println("\nModels with >1% error: $n_over1pct")
    println("\n--- Categories (not mutually exclusive) ---")
    @printf("  A ≠ U+F+C:           %d\n", length(cat_A_mismatch))
    @printf("  A = U+F+C but error:  %d\n", length(cat_A_match_err))
    @printf("  Negative entries:     %d\n", length(cat_negative))
    @printf("  R lambda ≈ 0:         %d\n", length(cat_zero_lambda))
    @printf("  R lambda > 100:       %d\n", length(cat_huge_lambda))
    @printf("  All-zero A:           %d\n", length(cat_all_zero))
    @printf("  Fixed by using A:     %d\n", length(cat_numerical))

    # Show some examples of A ≠ U+F+C
    if !isempty(cat_A_mismatch)
        println("\n--- Examples: A ≠ U+F+C (first 10) ---")
        for idx in cat_A_mismatch[1:min(10, length(cat_A_mismatch))]
            Ai = db.A[idx]; UFC = db.U[idx] .+ db.F[idx] .+ db.C[idx]
            diff = maximum(abs.(Ai .- UFC))
            r_lam = db.lambda_R[idx]
            jl_A = maximum(real.(eigen(Ai).values))
            jl_UFC = maximum(real.(eigen(UFC).values))
            @printf("  [%d] %s dim=%d  R=%.6f  Julia(A)=%.6f  Julia(U+F+C)=%.6f  max|A-(U+F+C)|=%.2e\n",
                    idx, db.species[idx], db.dims[idx], r_lam, jl_A, jl_UFC, diff)
        end
    end

    # Show examples where A=U+F+C but still error
    if !isempty(cat_A_match_err)
        println("\n--- Examples: A = U+F+C but >1% error (first 10) ---")
        for idx in cat_A_match_err[1:min(10, length(cat_A_match_err))]
            r_lam = db.lambda_R[idx]
            jl_lam = maximum(real.(eigen(db.A[idx]).values))
            rel = abs(jl_lam - r_lam) / max(abs(r_lam), 1e-10)
            @printf("  [%d] %s dim=%d  R=%.6f  Julia=%.6f  rel_err=%.2e\n",
                    idx, db.species[idx], db.dims[idx], r_lam, jl_lam, rel)
        end
    end

    # Overall: how many errors go away if we use stored A?
    valid_A = filter(isfinite, errors_from_A)
    valid_UFC = filter(isfinite, errors_from_UFC)
    println("\n--- Error from stored A vs from U+F+C ---")
    @printf("  From A:     n=%d  mean=%.2e  median=%.2e  >1%%=%d\n",
            length(valid_A), mean(valid_A), median(valid_A),
            count(x -> x >= 0.01, valid_A))
    @printf("  From U+F+C: n=%d  mean=%.2e  median=%.2e  >1%%=%d\n",
            length(valid_UFC), mean(valid_UFC), median(valid_UFC),
            count(x -> x >= 0.01, valid_UFC))

    return (cat_A_mismatch=cat_A_mismatch, cat_A_match_err=cat_A_match_err,
            cat_numerical=cat_numerical, n_over1pct=n_over1pct)
end

for db_name in ["compadre", "comadre"]
    try
        diagnose(db_name)
    catch e
        println("ERROR ($db_name): $e")
        println(sprint(showerror, e, catch_backtrace()))
    end
end
