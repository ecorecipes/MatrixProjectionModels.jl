@testset "Sampling Error" begin
    U = [0.0 0.0 0.0; 0.5 0.3 0.0; 0.0 0.4 0.2]
    F = [0.0 0.5 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    mpm = MatrixProjectionModel(U, F)
    rng = MersenneTwister(42)

    @testset "add_mpm_error" begin
        mpm_err = add_mpm_error(mpm, 100; rng=rng)
        @test n_stages(mpm_err) == 3
        @test all(mpm_err.A .>= 0)
        # With large sample, should be reasonably close to original
        mpm_err_large = add_mpm_error(mpm, 100000; rng=MersenneTwister(1))
        @test mpm_err_large.A ≈ mpm.A atol=0.15
    end

    @testset "add_mpm_error on raw matrix" begin
        A = U .+ F
        mpm_err = add_mpm_error(A, 100; rng=rng)
        @test n_stages(mpm_err) == 3
    end

    @testset "calculate_errors" begin
        C = [0.0 0.1 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        errs = calculate_errors(MatrixProjectionModel(U, F, C), 50; n_boot=100, rng=rng)
        @test size(errs.A) == (3, 3)
        @test size(errs.C) == (3, 3)
        @test all(errs.A .>= 0)
        @test errs.C[1, 2] ≥ 0
    end

    @testset "compute_ci" begin
        ci = compute_ci(mpm, m -> lambda(m), 50; n_boot=200, ci=0.99, rng=MersenneTwister(999))
        @test ci.lower < ci.upper
        @test ci.se > 0
        @test length(ci.boot_values) <= 200
    end
end
