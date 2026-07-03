@testset "Time Lag" begin

    # Standard Leslie matrix: A = U + F
    U = [0.0 0.0 0.0;
         0.5 0.0 0.0;
         0.0 0.3 0.0]
    F = [0.0 3.0 1.0;
         0.0 0.0 0.0;
         0.0 0.0 0.0]

    @testset "LaggedMPM construction" begin
        mpm = MatrixProjectionModel(U, F)
        lm = LaggedMPM(mpm; F_lag=1)
        @test lm isa LaggedMPM
        @test lm.lag_structure.max_lag == 1
        @test lm.U_lag == 0
        @test lm.F_lag == 1
        @test size(lm.augmented) == (6, 6)

        # Top-left block should be U (survival at lag 0)
        @test lm.augmented[1:3, 1:3] ≈ U
        # Top-right block should be F (fecundity at lag 1)
        @test lm.augmented[1:3, 4:6] ≈ F
        # Bottom-left block should be I
        @test lm.augmented[4:6, 1:3] ≈ Matrix{Float64}(I, 3, 3)
        # Bottom-right block should be 0
        @test lm.augmented[4:6, 4:6] ≈ zeros(3, 3)
    end

    @testset "LaggedMPM from kernel vector" begin
        K0 = [0.5 0.0; 0.3 0.2]
        K1 = [0.0 1.0; 0.0 0.0]
        lm = LaggedMPM([K0, K1])
        @test lm.lag_structure.max_lag == 1
        @test size(lm.augmented) == (4, 4)
        @test lm.augmented[1:2, 1:2] ≈ K0
        @test lm.augmented[1:2, 3:4] ≈ K1
    end

    @testset "LaggedMPM analysis" begin
        mpm = MatrixProjectionModel(U, F)
        lm = LaggedMPM(mpm; F_lag=1)

        λ = lambda(lm)
        @test λ > 0
        @test isfinite(λ)

        # Lagged λ should differ from standard λ
        λ_standard = lambda(mpm)
        @test !isapprox(λ, λ_standard; atol=1e-4)

        # Stable distribution should exist and be normalized
        w = stable_distribution(lm)
        @test length(w) == 6  # augmented dimension
        @test sum(w) ≈ 1.0 atol=1e-6

        # Damping ratio
        ρ = damping_ratio(lm)
        @test ρ > 1.0
    end

    @testset "LaggedMPM R0" begin
        mpm = MatrixProjectionModel(U, F)
        lm = LaggedMPM(mpm; F_lag=1)
        R0 = net_repro_rate(lm)
        @test R0 > 0
        @test isfinite(R0)
    end

    @testset "DirectIteration with LaggedMPM" begin
        mpm = MatrixProjectionModel(U, F)
        lm = LaggedMPM(mpm; F_lag=1)
        n0 = [10.0, 5.0, 2.0]
        prob = MPMProblem(lm, n0, (0, 200))
        sol = solve(prob, DirectIteration())

        @test sol.retcode == :Success
        @test length(sol.t) == 201
        @test length(sol.u) == 201
        # Physical state should have same dimension as n0
        @test length(sol.u[end]) == 3
        # Per-step lambdas should converge
        @test length(sol.lambdas) == 200

        # Converged lambda should match eigenanalysis
        λ_iter = sol.lambdas[end]
        λ_eigen = lambda(lm)
        @test λ_iter ≈ λ_eigen atol=0.01
    end

    @testset "EigenAnalysis with LaggedMPM" begin
        mpm = MatrixProjectionModel(U, F)
        lm = LaggedMPM(mpm; F_lag=1)
        n0 = [10.0, 5.0, 2.0]
        prob = MPMProblem(lm, n0, (0, 50))
        sol = solve(prob, EigenAnalysis())

        @test sol.retcode == :Success
        @test sol.eigenanalysis !== nothing
        @test sol.eigenanalysis.lambda > 0
    end

    @testset "Multi-lag L=2" begin
        mpm = MatrixProjectionModel(U, F)
        lm = LaggedMPM(mpm; U_lag=0, F_lag=2)
        @test lm.lag_structure.max_lag == 2
        @test size(lm.augmented) == (9, 9)

        λ = lambda(lm)
        @test λ > 0
        @test isfinite(λ)

        # DirectIteration should work
        n0 = [10.0, 5.0, 2.0]
        prob = MPMProblem(lm, n0, (0, 50))
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
    end

    @testset "F=0 degeneration" begin
        # With F_lag=1 but F=0, should behave like standard U-only model
        U_only = [0.0 0.0; 0.5 0.0]
        mpm = MatrixProjectionModel(U_only)
        lm = LaggedMPM(mpm; F_lag=1)
        @test lambda(lm) ≈ lambda(U_only) atol=1e-10
    end

    @testset "Normalized iteration" begin
        mpm = MatrixProjectionModel(U, F)
        lm = LaggedMPM(mpm; F_lag=1)
        n0 = [10.0, 5.0, 2.0]
        prob = MPMProblem(LeslieMPM(), DensityIndependent(), Deterministic(),
                          lm, n0, (0, 20); normalize=true)
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        # With normalization, population should stay bounded
        @test all(sum.(sol.u) .< 1e10)
    end

    @testset "Extinction guard" begin
        lm = LaggedMPM([zeros(2, 2), zeros(2, 2)])
        sol = solve(MPMProblem(lm, [1.0, 0.0], (0, 3)), DirectIteration())
        @test all(isfinite, sol.lambdas)
    end

end
