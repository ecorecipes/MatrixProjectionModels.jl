# Life History Analysis
# Demonstrates Rage-equivalent analysis of a stage-classified MPM

using MatrixProjectionModels

# A 3-stage MPM: juvenile → adult → senescent
U = [0.0 0.0 0.0;
     0.5 0.3 0.0;
     0.0 0.4 0.2]
F = [0.0 0.5 3.0;
     0.0 0.0 0.0;
     0.0 0.0 0.0]

mpm = MatrixProjectionModel(U, F)
A = mpm.A

println("=== Matrix Population Model ===")
println("λ = ", round(lambda(A); digits=4))

# Life table
println("\n=== Life Table ===")
tbl = mpm_to_table(U, F)
for i in 1:min(10, length(tbl.x))
    println("  Age $(tbl.x[i]): lx=$(round(tbl.lx[i]; digits=3)) px=$(round(tbl.px[i]; digits=3)) mx=$(round(tbl.mx[i]; digits=3))")
end

# Life history traits
println("\n=== Life History Traits ===")
println("Mean life expectancy: ", round(life_expect_mean(U); digits=2))
println("Longevity (lx < 0.01): ", longevity(U))
println("Prob of maturity: ", round(mature_prob(U, F); digits=3))
println("Mean age at maturity: ", round(mature_age(U, F); digits=2))
println("R₀ (net repro rate): ", round(net_repro_rate(U, F); digits=3))
println("Generation time (T): ", round(gen_time(U, F); digits=2))

# Entropy
lx = mpm_to_lx(U)
println("Keyfitz's entropy: ", round(entropy_k(lx); digits=3))
println("Survivorship shape: ", round(shape_surv(lx); digits=4))

# Vital rates
println("\n=== Vital Rates ===")
println("Survival: ", round.(vr_vec_survival(U); digits=3))
println("Growth: ", round.(vr_vec_growth(U); digits=3))
println("Stasis: ", round.(vr_vec_stasis(U); digits=3))
println("Mean survival: ", round(vr_survival(U); digits=3))

# Perturbation
println("\n=== Perturbation ===")
E = elasticity(A)
println("Elasticity matrix:")
display(round.(E; digits=3))

vr_pert = perturb_vr(U, F; type=:elasticity)
println("\nVital rate elasticities:")
for (k, v) in pairs(vr_pert)
    println("  $k: $(round(v; digits=4))")
end
