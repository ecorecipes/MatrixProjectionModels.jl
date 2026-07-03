@testset "Density Response" begin
    A = [1.0 2.0; 3.0 4.0]
    spec = DensityVitalRateSpec(
        survival = BevertonHoltDensity(α=2.0, β=0.1),
        fecundity = RickerDensity(α=0.0, β=0.2),
    )
    N = 5.0

    @testset "combined matrix default fecundity rows" begin
        A_dd = apply_density(spec, A, N)
        @test A_dd[1, :] ≈ A[1, :] .* spec.fecundity(N)
        @test A_dd[2, :] ≈ A[2, :] .* spec.survival(N)
    end

    @testset "unsupported density channels fail loudly" begin
        @test_throws ArgumentError apply_density(
            DensityVitalRateSpec(growth=LogisticDensity(K=10.0)),
            A, N
        )
        @test_throws ArgumentError apply_density(
            DensityVitalRateSpec(recruitment=UsherDensity(α=0.0, β=0.1)),
            A, N
        )
        @test_throws ArgumentError apply_density(
            DensityVitalRateSpec(time_delay=2),
            A, N
        )
    end
end
