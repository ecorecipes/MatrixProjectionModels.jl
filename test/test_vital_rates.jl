@testset "Vital Rates" begin
    # Lefkovitch-style matrix with stasis, growth, and retrogression
    U = [0.2 0.0 0.0;
         0.3 0.4 0.1;
         0.0 0.2 0.3]
    F = [0.0 1.0 3.0;
         0.0 0.0 0.0;
         0.0 0.0 0.0]

    @testset "vr_vec_survival" begin
        surv = vr_vec_survival(U)
        @test length(surv) == 3
        @test surv[1] ≈ 0.5  # 0.2 + 0.3
        @test surv[2] ≈ 0.6  # 0.4 + 0.2
        @test surv[3] ≈ 0.4  # 0.1 + 0.3
    end

    @testset "vr_vec_growth" begin
        growth = vr_vec_growth(U)
        @test length(growth) == 3
        # Stage 1: growth = transition to stage 1 (higher = lower index)
        # Growth means moving to a "higher" (lower index) stage
        @test growth[1] ≈ 0.0  # Can't grow from stage 1
        @test growth[2] ≈ 0.0 / 0.6  # No transition from 2 to 1 through upper triangle
    end

    @testset "vr_vec_shrinkage" begin
        shrinkage = vr_vec_shrinkage(U)
        @test length(shrinkage) == 3
        @test shrinkage[1] ≈ 0.3 / 0.5  # U[2,1]/sum = transition to lower stage
    end

    @testset "vr_vec_stasis" begin
        stasis = vr_vec_stasis(U)
        @test stasis[1] ≈ 0.2 / 0.5
        @test stasis[2] ≈ 0.4 / 0.6
        @test stasis[3] ≈ 0.3 / 0.4
    end

    @testset "vr_vec_reproduction" begin
        repro = vr_vec_reproduction(U, F)
        @test repro[1] ≈ 0.0  # No fecundity from stage 1
        @test repro[2] ≈ 1.0 / 0.6  # F column sum / survival
        @test repro[3] ≈ 3.0 / 0.4
    end

    @testset "Weighted averages" begin
        s = vr_survival(U)
        @test 0 < s < 1

        g = vr_growth(U)
        @test g >= 0

        f = vr_fecundity(U, F)
        @test f >= 0
    end

    @testset "Dormancy" begin
        # Mark stage 3 as dormant
        dorm = [3]
        enter = vr_vec_dorm_enter(U; dorm_stages=dorm)
        exit_d = vr_vec_dorm_exit(U; dorm_stages=dorm)
        @test enter[3] == 0.0  # Stage 3 is dormant, no "entering" from itself
        @test exit_d[1] == 0.0  # Stage 1 is not dormant
    end
end
