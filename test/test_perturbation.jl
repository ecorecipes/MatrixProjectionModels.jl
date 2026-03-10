@testset "Perturbation" begin
    U = [0.0 0.0 0.0; 0.5 0.3 0.0; 0.0 0.4 0.2]
    F = [0.0 0.5 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    A = U .+ F

    @testset "perturb_matrix sensitivity" begin
        S = perturb_matrix(A; type=:sensitivity)
        @test size(S) == (3, 3)
        @test all(isfinite.(S))
        # Should be close to analytical sensitivity
        S_analytical = sensitivity(A)
        @test S ≈ S_analytical atol=1e-3
    end

    @testset "perturb_matrix elasticity" begin
        E = perturb_matrix(A; type=:elasticity)
        @test size(E) == (3, 3)
        # Non-zero elements elasticity should sum ≈ 1
        e_sum = sum(E[i, j] for i in 1:3, j in 1:3 if A[i, j] != 0)
        @test e_sum ≈ 1.0 atol=0.01
    end

    @testset "perturb_matrix custom stat" begin
        S = perturb_matrix(A; type=:sensitivity, demog_stat=A -> damping_ratio(A))
        @test all(isfinite.(S))
    end

    @testset "perturb_vr" begin
        result = perturb_vr(U, F; type=:sensitivity)
        @test haskey(result, :survival)
        @test haskey(result, :growth)
        @test haskey(result, :fecundity)
        @test all(isfinite(v) for v in values(result))
    end

    @testset "perturb_trans" begin
        result = perturb_trans(U, F; type=:sensitivity)
        @test size(result.U) == (3, 3)
        @test size(result.F) == (3, 3)
        @test size(result.C) == (3, 3)
    end
end
