@testset "Analysis" begin
    # Classic 2-stage Leslie matrix
    A = [0.0 3.0; 0.5 0.0]
    λ_exact = sqrt(1.5)  # ≈ 1.2247

    @testset "lambda" begin
        @test lambda(A) ≈ λ_exact atol=1e-6
        mpm = MatrixProjectionModel(A)
        @test lambda(mpm) ≈ λ_exact atol=1e-6
    end

    @testset "stable_distribution" begin
        w = stable_distribution(A)
        @test length(w) == 2
        @test sum(w) ≈ 1.0 atol=1e-10
        @test all(w .>= 0)
        # Should satisfy A*w = λ*w (up to normalization)
        Aw = A * w
        @test Aw / sum(Aw) ≈ w atol=1e-6
    end

    @testset "reproductive_value" begin
        v = reproductive_value(A)
        @test length(v) == 2
        # dot(v, w) should be 1
        w = stable_distribution(A)
        @test dot(v, w) ≈ 1.0 atol=1e-6
    end

    @testset "sensitivity" begin
        S = sensitivity(A)
        @test size(S) == (2, 2)
        @test all(isfinite.(S))
        # Sensitivity of λ to A[i,j] = v[i]*w[j]/dot(v,w)
        # All sensitivities should be positive for this matrix
        @test all(S .> 0)
    end

    @testset "elasticity" begin
        E = elasticity(A)
        @test size(E) == (2, 2)
        # Elasticities of non-zero elements should sum to 1
        e_sum = sum(E[i, j] for i in 1:2, j in 1:2 if A[i, j] != 0)
        @test e_sum ≈ 1.0 atol=1e-6
    end

    @testset "damping_ratio" begin
        dr = damping_ratio(A)
        @test dr > 1.0  # Dominant > subdominant
        @test isfinite(dr)
    end

    @testset "3-stage matrix" begin
        # Rage's mpm1 test data equivalent
        U = [0.0 0.0 0.0; 0.5 0.3 0.0; 0.0 0.4 0.2]
        F = [0.0 0.5 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        A = U .+ F
        @test lambda(A) > 0
        w = stable_distribution(A)
        @test sum(w) ≈ 1.0 atol=1e-10
        E = elasticity(A)
        e_sum = sum(E[i, j] for i in 1:3, j in 1:3 if A[i, j] != 0)
        @test e_sum ≈ 1.0 atol=1e-6
    end
end
