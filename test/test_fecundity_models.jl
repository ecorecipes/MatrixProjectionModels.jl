@testset "Fecundity Models" begin
    @testset "LogisticFecundity" begin
        m = LogisticFecundity(5.0, 0.5, 10.0)
        @test m(10.0) ≈ 2.5  # Midpoint
        @test m(0.0) < 0.5   # Low at young age
        @test m(100.0) ≈ 5.0 atol=0.01  # Asymptote
    end

    @testset "StepFecundity" begin
        m = StepFecundity(3.0, 5.0)
        @test m(4.0) == 0.0
        @test m(5.0) == 3.0
        @test m(10.0) == 3.0
    end

    @testset "VonBertalanffyFecundity" begin
        m = VonBertalanffyFecundity(10.0, 0.1)
        @test m(0.0) ≈ 0.0
        @test m(100.0) ≈ 10.0 atol=0.1
        @test m(10.0) > 0
    end

    @testset "NormalFecundity" begin
        m = NormalFecundity(5.0, 20.0, 3.0)
        @test m(20.0) ≈ 5.0  # Peak at mu
        @test m(0.0) < 0.01  # Far from peak
    end

    @testset "HadwigerFecundity" begin
        m = HadwigerFecundity(1.0, 3.5, 25.0)
        @test m(0.0) ≈ 0.0
        @test m(25.0) > 0  # Peak around C
    end

    @testset "model_fecundity" begin
        m = LogisticFecundity(5.0, 0.5, 10.0)
        result = model_fecundity(m; ages=0:50)
        @test length(result.x) == 51
        @test length(result.fx) == 51
        @test result.fx[11] ≈ 2.5 atol=0.1
    end
end
