@testset "Transformation" begin
    U = [0.0 0.0 0.0; 0.5 0.3 0.0; 0.0 0.4 0.2]
    F = [0.0 0.5 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0]

    @testset "repro_stages" begin
        r = repro_stages(F)
        @test r == [false, true, true]
    end

    @testset "standard_stages" begin
        stages = standard_stages(F)
        @test length(stages) == 3
        @test stages[1] == PropaguleStage  # Receives reproduction but doesn't reproduce
    end

    @testset "mpm_split" begin
        A = U .+ F
        mpm = mpm_split(A)
        @test n_stages(mpm) == 3
        # First row should be in F
        @test mpm.F[1, :] ≈ A[1, :]
        # Rest in U
        @test mpm.U[1, :] ≈ zeros(3)
        @test mpm.U[2:end, :] ≈ A[2:end, :]
    end

    @testset "mpm_collapse" begin
        # Collapse stages 1+2 into one, keep 3
        collapse = [[1, 2], [3]]
        mpm_c = mpm_collapse(U, F, collapse)
        @test n_stages(mpm_c) == 2

        # Lambda should be approximately preserved
        A = U .+ F
        @test lambda(mpm_c) ≈ lambda(A) atol=0.5
    end

    @testset "mpm_standardize" begin
        mpm_std = mpm_standardize(U, F)
        @test n_stages(mpm_std) == 3
        # All reproduction should be in first row
        @test all(mpm_std.F[2:end, :] .== 0)
        @test sum(mpm_std.F[1, :]) ≈ sum(F) atol=1e-10
    end

    @testset "mpm_rearrange" begin
        mpm_r = mpm_rearrange(U, F; new_order=[3, 2, 1])
        @test n_stages(mpm_r) == 3
        # Rearranging should permute rows and columns
        @test mpm_r.U[1, 1] ≈ U[3, 3]
        @test mpm_r.F[3, 3] ≈ F[1, 1]
    end

    @testset "mpm_rearrange validation" begin
        @test_throws ArgumentError mpm_rearrange(U, F; new_order=[1, 2])
        @test_throws ArgumentError mpm_rearrange(U, F; new_order=[1, 1, 2])
    end

    @testset "name_stages" begin
        mpm = MatrixProjectionModel(U, F)
        mpm_named = name_stages(mpm, ["juvenile", "adult", "senescent"])
        @test mpm_named.stage_names == [:juvenile, :adult, :senescent]
    end
end
