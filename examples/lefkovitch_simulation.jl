# Lefkovitch Matrix Simulation
# Demonstrates stochastic simulation with stage-structured MPMs

using MatrixProjectionModels
using Random

rng = MersenneTwister(42)

# 1. Generate random Lefkovitch matrices (different archetypes)
println("Archetype comparison:")
for arch in 1:4
    mpm = rand_lefko_mpm(4, [0.0, 0.5, 2.0, 3.0]; archetype=arch, rng=MersenneTwister(42))
    println("  Archetype $arch: λ = $(round(lambda(mpm); digits=4))")
end

# 2. Deterministic simulation
mpm = rand_lefko_mpm(4, [0.0, 0.5, 2.0, 3.0]; archetype=3, rng=rng)
n0 = [100.0, 50.0, 30.0, 10.0]
prob = MPMProblem(mpm, n0, (0, 50))
sol = solve(prob, DirectIteration())

println("\nDeterministic simulation (50 steps):")
println("  Initial N: $(sum(n0))")
println("  Final N: $(round(sum(sol.u[end]); digits=1))")
println("  λ (eigenanalysis): $(round(lambda(sol); digits=6))")

# 3. Stochastic simulation with kernel resampling
mpms = rand_lefko_set(5; n_stages=4, fecundity=[0.0, 0.5, 2.0, 3.0],
                      archetype=3, rng=rng)
prob_stoch = MPMProblem(mpms, n0, (0, 1000))
sol_stoch = solve(prob_stoch, DirectIteration(); rng=rng)

log_λs = stochastic_growth_rate(sol_stoch; burn_in=100)
println("\nStochastic simulation (1000 steps, 5 environments):")
println("  Stochastic λ: $(round(log_λs; digits=6))")
