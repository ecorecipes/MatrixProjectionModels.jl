using Distributions

@testset "Leslie Construction" begin
    @testset "make_leslie_mpm from vectors" begin
        surv = [0.8, 0.5, 0.3]
        fec = [0.0, 1.0, 2.5, 1.5]
        mpm = make_leslie_mpm(surv, fec)

        @test n_stages(mpm) == 4
        @test is_leslie(mpm.A)
        @test mpm.F[1, :] ≈ fec
        @test mpm.U[2, 1] ≈ 0.8
        @test mpm.U[3, 2] ≈ 0.5
        @test mpm.U[4, 3] ≈ 0.3
        @test mpm.A ≈ mpm.U .+ mpm.F
    end

    @testset "make_leslie_mpm from models" begin
        mort = GompertzMortality(0.01, 0.1)
        fec = LogisticFecundity(5.0, 0.5, 10.0)
        mpm = make_leslie_mpm(mort, fec; truncate=0.01)

        @test n_stages(mpm) >= 2
        @test is_leslie(mpm.A)
        @test lambda(mpm) > 0
    end

    @testset "make_leslie_mpm with n_stages" begin
        mort = ExponentialMortality(0.1)
        fec = StepFecundity(2.0, 3.0)
        mpm = make_leslie_mpm(mort, fec; n_stages=10)

        @test n_stages(mpm) == 10
    end

    @testset "dimension validation" begin
        @test_throws ArgumentError make_leslie_mpm([0.5], [1.0, 2.0, 3.0])
    end

    @testset "rand_leslie_set samples parameters reproducibly" begin
        mort = ExponentialMortality(0.1)
        fec = LogisticFecundity(4.0, 0.4, 8.0)
        mort_dist = (C = Uniform(0.05, 0.15),)
        fec_dist = (A = Uniform(3.0, 5.0), k = Uniform(0.2, 0.6), x_mid = Uniform(7.0, 9.0))

        set1 = rand_leslie_set(4;
            mortality_model = mort,
            fecundity_model = fec,
            mortality_params_dist = mort_dist,
            fecundity_params_dist = fec_dist,
            n_stages = 8,
            rng = MersenneTwister(123))
        set2 = rand_leslie_set(4;
            mortality_model = mort,
            fecundity_model = fec,
            mortality_params_dist = mort_dist,
            fecundity_params_dist = fec_dist,
            n_stages = 8,
            rng = MersenneTwister(123))

        @test any(!isapprox(set1[1].A, m.A) for m in set1[2:end])
        @test all(set1[i].A ≈ set2[i].A for i in eachindex(set1, set2))
    end
end
