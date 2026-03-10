"""
Compare R and Julia MPM benchmark results.
Run after benchmark_r.R and benchmark_julia.jl.
"""

using Printf

# R results (manually extracted from benchmark output)
r_compadre = (
    lambda = 0.0293,
    sensitivity = 0.0593,
    elasticity = 0.0967,
    life_expect = 0.0077,
    net_repro_rate = 0.0087,
    gen_time = 0.0267,
    longevity = 0.0167,
    full_pipeline = 0.2323,
)

r_comadre = (
    lambda = 0.0290,
    sensitivity = 0.0530,
    elasticity = 0.0810,
    life_expect = 0.0077,
    net_repro_rate = 0.0070,
    gen_time = 0.0277,
    longevity = 0.0137,
    full_pipeline = 0.2143,
)

# Julia results (after longevity optimization)
jl_compadre = (
    lambda = 0.0039,
    sensitivity = 0.0080,
    elasticity = 0.0078,
    life_expect = 0.0013,
    net_repro_rate = 0.0013,
    gen_time = 0.0123,
    longevity = 0.0021,
    full_pipeline = 0.0461,
)

jl_comadre = (
    lambda = 0.0023,
    sensitivity = 0.0050,
    elasticity = 0.0049,
    life_expect = 0.0017,
    net_repro_rate = 0.0005,
    gen_time = 0.0095,
    longevity = 0.0015,
    full_pipeline = 0.0391,
)

# Single matrix (R amortized, Julia BenchmarkTools median)
r_single = Dict(
    "lambda 3x3" => 46.3, "sensitivity 3x3" => 99.8, "full 3x3" => 394.3,
    "lambda 5x5" => 48.9, "sensitivity 5x5" => 103.2, "full 5x5" => 424.0,
    "lambda 10x10" => 58.3, "sensitivity 10x10" => 120.6, "full 10x10" => 494.8,
)
jl_single = Dict(
    "lambda 3x3" => 2.5, "sensitivity 3x3" => 5.4, "full 3x3" => 25.5,
    "lambda 5x5" => 5.4, "sensitivity 5x5" => 11.3, "full 5x5" => 49.5,
    "lambda 10x10" => 12.8, "sensitivity 10x10" => 26.2, "full 10x10" => 113.3,
)

println("=" ^ 75)
println("  Julia (MatrixProjectionModels.jl) vs R — Benchmark Comparison")
println("  500 models from COMPADRE, same hardware, same data")
println("=" ^ 75)

println("\n--- Batch benchmarks: 500 COMPADRE models ---")
println("-" ^ 75)
@printf("%-20s  %10s  %10s  %10s\n", "Function", "R (s)", "Julia (s)", "Speedup")
println("-" ^ 75)

funcs = [:lambda, :sensitivity, :elasticity, :life_expect,
         :net_repro_rate, :gen_time, :longevity, :full_pipeline]
for f in funcs
    r = getfield(r_compadre, f)
    j = getfield(jl_compadre, f)
    speedup = r / j
    marker = speedup >= 1.0 ? "" : " (R faster)"
    @printf("%-20s  %10.4f  %10.4f  %8.1fx%s\n", f, r, j, speedup, marker)
end

println("\n--- Batch benchmarks: 500 COMADRE models ---")
println("-" ^ 75)
@printf("%-20s  %10s  %10s  %10s\n", "Function", "R (s)", "Julia (s)", "Speedup")
println("-" ^ 75)
for f in funcs
    r = getfield(r_comadre, f)
    j = getfield(jl_comadre, f)
    speedup = r / j
    marker = speedup >= 1.0 ? "" : " (R faster)"
    @printf("%-20s  %10.4f  %10.4f  %8.1fx%s\n", f, r, j, speedup, marker)
end

println("\n--- Single matrix benchmarks (us/call) ---")
println("-" ^ 75)
@printf("%-22s  %10s  %10s  %10s\n", "Operation", "R (us)", "Julia (us)", "Speedup")
println("-" ^ 75)
for key in ["lambda 3x3", "sensitivity 3x3", "full 3x3",
            "lambda 5x5", "sensitivity 5x5", "full 5x5",
            "lambda 10x10", "sensitivity 10x10", "full 10x10"]
    r = r_single[key]
    j = jl_single[key]
    speedup = r / j
    @printf("%-22s  %10.1f  %10.1f  %8.1fx\n", key, r, j, speedup)
end
println("-" ^ 75)
