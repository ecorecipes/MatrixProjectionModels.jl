@testset "Mortality Models" begin
    @testset "GompertzMortality" begin
        m = GompertzMortality(0.01, 0.1)
        @test m(0) ≈ 0.01
        @test m(10) ≈ 0.01 * exp(1.0)
        @test m(20) > m(10)  # Increasing hazard
    end

    @testset "GompertzMakehamMortality" begin
        m = GompertzMakehamMortality(0.01, 0.1, 0.005)
        @test m(0) ≈ 0.01 + 0.005
        @test m(10) ≈ 0.01 * exp(1.0) + 0.005
    end

    @testset "ExponentialMortality" begin
        m = ExponentialMortality(0.1)
        @test m(0) ≈ 0.1
        @test m(100) ≈ 0.1  # Constant
    end

    @testset "SilerMortality" begin
        m = SilerMortality(0.1, 0.5, 0.01, 0.001, 0.1)
        @test m(0) ≈ 0.1 + 0.01 + 0.001
        # Bathtub shape: high at 0, lower in middle, high at old age
        h0 = m(0)
        h10 = m(10)
        h50 = m(50)
        @test h10 < h0  # Juvenile component decreases
        @test h50 > h10  # Senescent component increases
    end

    @testset "WeibullMortality" begin
        m = WeibullMortality(2.0, 0.1)
        @test m(0) ≈ 0.0  # h(0) = 0 for b0 > 1
        @test m(10) > 0
    end

    @testset "WeibullMakehamMortality" begin
        m = WeibullMakehamMortality(2.0, 0.1, 0.05)
        @test m(0) ≈ 0.05  # Only Makeham component at x=0
    end

    @testset "model_survival" begin
        m = GompertzMortality(0.01, 0.1)
        result = model_survival(m; truncate=0.01)
        @test result.lx[1] ≈ 1.0 atol=0.1  # Approximately 1 at age 0
        @test result.lx[end] < 0.01 + 0.01  # Truncated
        @test length(result.x) == length(result.hx) == length(result.lx)
        @test length(result.px) == length(result.lx)
        @test all(result.hx .>= 0)
        @test all(0 .<= result.px .<= 1)
    end
end
