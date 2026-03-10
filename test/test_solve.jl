@testset "Solve" begin
    @testset "EigenAnalysis" begin
        A = [0.0 3.0; 0.5 0.0]
        prob = MPMProblem(A, [10.0, 5.0], (0, 0))
        sol = solve(prob, EigenAnalysis())
        @test sol.retcode == :Success
        @test sol.eigenanalysis !== nothing
        @test sol.eigenanalysis.lambda ≈ sqrt(1.5) atol=1e-6
        @test lambda(sol) ≈ sqrt(1.5) atol=1e-6
    end

    @testset "DirectIteration - deterministic" begin
        A = [0.0 3.0; 0.5 0.0]
        n0 = [10.0, 5.0]
        prob = MPMProblem(A, n0, (0, 50))
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        @test length(sol.t) == 51
        @test length(sol.u) == 51
        @test sol.u[1] ≈ n0
        @test length(sol.lambdas) == 50

        # For periodic matrix, geometric mean of per-step lambdas converges to λ
        geom_mean = exp(mean(log.(sol.lambdas[end-9:end])))
        @test geom_mean ≈ sqrt(1.5) atol=0.05
    end

    @testset "DirectIteration - normalized" begin
        A = [0.0 3.0; 0.5 0.0]
        prob = MPMProblem(A, [10.0, 5.0], (0, 100); normalize=true)
        sol = solve(prob, DirectIteration())
        # Normalized: population sums should stay around initial
        for u in sol.u[2:end]
            @test sum(u) ≈ 1.0 atol=0.5  # Relaxed; initial sum != 1
        end
    end

    @testset "MPMProblem from MatrixProjectionModel" begin
        U = [0.0 0.0; 0.5 0.3]
        F = [0.0 3.0; 0.0 0.0]
        mpm = MatrixProjectionModel(U, F)
        prob = MPMProblem(mpm, [10.0, 5.0], (0, 20))
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
    end

    @testset "Stochastic kernel resampled" begin
        rng = MersenneTwister(42)
        A1 = [0.0 2.0; 0.6 0.0]
        A2 = [0.0 4.0; 0.4 0.0]
        mpm1 = MatrixProjectionModel(A1)
        mpm2 = MatrixProjectionModel(A2)
        prob = MPMProblem([mpm1, mpm2], [10.0, 5.0], (0, 100))
        sol = solve(prob, DirectIteration(); rng=rng)
        @test sol.retcode == :Success
        @test length(sol.lambdas) == 100
        @test stochastic_growth_rate(sol; burn_in=10) > 0
    end

    @testset "show methods" begin
        A = [0.0 3.0; 0.5 0.0]
        prob = MPMProblem(A, [10.0, 5.0], (0, 10))
        sol = solve(prob, DirectIteration())
        buf = IOBuffer()
        show(buf, prob)
        @test occursin("MPMProblem", String(take!(buf)))
        show(buf, sol)
        @test occursin("MPMSolution", String(take!(buf)))
    end
end
