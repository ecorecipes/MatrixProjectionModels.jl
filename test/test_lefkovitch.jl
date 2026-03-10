@testset "Lefkovitch Construction" begin
    rng = MersenneTwister(123)

    @testset "Archetype $arch" for arch in 1:4
        mpm = rand_lefko_mpm(4, [0.0, 0.5, 2.0, 3.0]; archetype=arch, rng=MersenneTwister(42))
        @test n_stages(mpm) == 4
        @test mpm.F[1, :] ≈ [0.0, 0.5, 2.0, 3.0]
        # Column sums of U should be <= 1 (survival probabilities)
        for j in 1:4
            @test sum(mpm.U[:, j]) <= 1.0 + 1e-10
            @test sum(mpm.U[:, j]) >= 0.0
        end
        @test lambda(mpm) > 0
    end

    @testset "Scalar fecundity" begin
        mpm = rand_lefko_mpm(3, 2.0; archetype=1, rng=rng)
        @test mpm.F[1, :] ≈ [2.0, 2.0, 2.0]
    end

    @testset "Validation" begin
        @test_throws ArgumentError rand_lefko_mpm(3, [1.0, 2.0]; archetype=1)
        @test_throws ArgumentError rand_lefko_mpm(3, 2.0; archetype=5)
        @test_throws ArgumentError rand_lefko_mpm(1, [2.0]; archetype=1)
    end

    @testset "rand_lefko_set" begin
        mpms = rand_lefko_set(5; n_stages=3, fecundity=2.0, archetype=2, rng=rng)
        @test length(mpms) == 5
        @test all(m -> n_stages(m) == 3, mpms)
    end

    @testset "rand_lefko_set with constraint" begin
        mpms = rand_lefko_set(3; n_stages=3, fecundity=2.0, archetype=1,
                              constraint=m -> lambda(m) > 1.0, rng=rng)
        @test all(m -> lambda(m) > 1.0, mpms)
    end
end
