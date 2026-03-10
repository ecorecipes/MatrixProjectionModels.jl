# Leslie Matrix Basics
# Demonstrates construction and analysis of a simple Leslie MPM

using MatrixProjectionModels

# 1. Build a Leslie matrix from survival and fecundity schedules
survival = [0.8, 0.7, 0.5, 0.3]
fecundity = [0.0, 0.5, 2.0, 3.0, 1.0]
mpm = make_leslie_mpm(survival, fecundity)

println("Projection matrix A:")
display(mpm.A)

# 2. Eigenanalysis
println("\nDominant eigenvalue (λ): ", lambda(mpm))
println("Stable stage distribution: ", stable_distribution(mpm))
println("Reproductive value: ", reproductive_value(mpm))
println("Damping ratio: ", damping_ratio(mpm))

# 3. Sensitivity and elasticity
println("\nSensitivity matrix:")
display(sensitivity(mpm))

println("\nElasticity matrix:")
E = elasticity(mpm)
display(E)
println("Sum of elasticities: ", sum(E[i,j] for i in 1:5, j in 1:5 if mpm.A[i,j] != 0))

# 4. Build from mortality/fecundity models
mort = GompertzMortality(0.01, 0.08)
fec = LogisticFecundity(5.0, 0.3, 15.0)

surv_schedule = model_survival(mort; truncate=0.01)
println("\nGompertz survival: $(length(surv_schedule.x)) age classes")

mpm2 = make_leslie_mpm(mort, fec)
println("λ from Gompertz/Logistic model: ", lambda(mpm2))
