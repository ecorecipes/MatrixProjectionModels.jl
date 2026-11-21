import ResilienceCore
using StructuredPopulationCore: DensityDependent

@testset "linearise: resilience operators" begin
    A = [0.0 3.0 1.0; 0.5 0.0 0.0; 0.0 0.3 0.0]
    m = MatrixProjectionModel(A)
    op = linearise(m)
    @test op isa ResilienceCore.ProjectionMatrix
    @test ResilienceCore.growth_rate(op) ≈ lambda(A)
    prob = MPMProblem(m, ones(3), (0, 10))
    @test ResilienceCore.growth_rate(linearise(prob)) ≈ lambda(A)
    @test ResilienceCore.resilience(linearise(prob)) ≈ log(damping_ratio(A))
    @test ResilienceCore.growth_rate(linearise(MPMProblem(A, ones(3), (0, 10)))) ≈ lambda(A)
    # no single projection matrix: density dependence, stochasticity
    @test_throws ArgumentError linearise(MPMProblem(DensityDependent(), (n, t, p) -> A, ones(3), (0, 10)))
    @test_throws ArgumentError linearise(MPMProblem([m, MatrixProjectionModel(0.9A)], ones(3), (0, 10)))
end
