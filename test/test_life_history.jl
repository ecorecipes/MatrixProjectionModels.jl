@testset "Life History" begin
    U = [0.0 0.0 0.0; 0.5 0.3 0.0; 0.0 0.4 0.2]
    F = [0.0 0.5 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    A = U .+ F

    @testset "life_expect_mean" begin
        le = life_expect_mean(U)
        @test le > 0
        @test isfinite(le)
        # Different starting stages
        le2 = life_expect_mean(U; start=2)
        @test le2 > 0
    end

    @testset "life_expect_var" begin
        v = life_expect_var(U)
        @test v > 0 || v ≈ 0  # Variance is non-negative
        @test isfinite(v)
    end

    @testset "mature_prob" begin
        prob = mature_prob(U, F)
        @test 0 <= prob <= 1
    end

    @testset "mature_age" begin
        age = mature_age(U, F)
        @test age >= 0
        @test isfinite(age)
    end

    @testset "net_repro_rate" begin
        R0 = net_repro_rate(U, F; method=:generation)
        @test R0 >= 0
        @test isfinite(R0)

        R0_fund = net_repro_rate(U, F; method=:fundamental)
        @test R0_fund >= 0
        # Both methods should give similar results
        @test R0 ≈ R0_fund atol=0.5
    end

    @testset "gen_time" begin
        T_R0 = gen_time(U, F; method=:R0)
        @test T_R0 > 0
        @test isfinite(T_R0)

        T_cohort = gen_time(U, F; method=:cohort)
        @test T_cohort > 0
        @test isfinite(T_cohort)

        T_age_diff = gen_time(U, F; method=:age_diff)
        @test T_age_diff > 0
        @test isfinite(T_age_diff)
    end

    @testset "longevity" begin
        l = longevity(U)
        @test l > 0
        @test l <= 1001  # Bounded by xmax
    end

    @testset "entropy_k" begin
        lx = mpm_to_lx(U)
        H = entropy_k(lx)
        @test H >= 0
        @test isfinite(H)

        H_age = entropy_k_age(U)
        @test H_age >= 0
    end

    @testset "entropy_d" begin
        lx = mpm_to_lx(U)
        mx = mpm_to_mx(U, F)
        s = entropy_d(lx, mx)
        @test s >= 0
        @test isfinite(s)
    end

    @testset "shape_surv" begin
        lx = mpm_to_lx(U)
        s = shape_surv(lx)
        @test -0.5 <= s <= 0.5
    end

    @testset "shape_rep" begin
        s = shape_rep(U, F)
        @test -0.5 <= s <= 0.5
    end

    @testset "mature_distrib" begin
        r_stages = repro_stages(F)
        dist = mature_distrib(U; repro_stages=r_stages)
        @test length(dist) == 3
        @test sum(dist) ≈ 1.0 atol=1e-6
    end
end
