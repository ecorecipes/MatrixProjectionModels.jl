@testset "Types" begin
    @testset "MatrixProjectionModel constructors" begin
        # From U + F
        U = [0.0 0.0; 0.5 0.3]
        F = [0.0 3.0; 0.0 0.0]
        mpm = MatrixProjectionModel(U, F)
        @test size(mpm) == (2, 2)
        @test mpm.A ≈ [0.0 3.0; 0.5 0.3]
        @test mpm.U ≈ U
        @test mpm.F ≈ F
        @test mpm.C ≈ zeros(2, 2)
        @test n_stages(mpm) == 2

        # From U + F + C
        C = [0.0 0.5; 0.0 0.0]
        mpm2 = MatrixProjectionModel(U, F, C)
        @test mpm2.A ≈ U .+ F .+ C
        @test mpm2.C ≈ C

        # From A only
        A = [0.0 3.0; 0.5 0.3]
        mpm3 = MatrixProjectionModel(A)
        @test mpm3.A ≈ A
        @test mpm3.U ≈ A  # When only A given, U = A

        # From survival/fecundity vectors (Leslie)
        surv = [0.8, 0.5]
        fec = [0.0, 1.5, 2.0]
        mpm4 = MatrixProjectionModel(surv, fec)
        @test mpm4.A[1, :] ≈ fec
        @test mpm4.A[2, 1] ≈ 0.8
        @test mpm4.A[3, 2] ≈ 0.5
        @test mpm4.A[2, 2] == 0.0
        @test mpm4.A[3, 3] == 0.0

        # Type promotion
        U_int = [0 0; 1 0]
        F_float = [0.0 3.0; 0.0 0.0]
        mpm5 = MatrixProjectionModel(U_int, F_float)
        @test eltype(mpm5) == Float64
    end

    @testset "AbstractMatrix interface" begin
        A = [0.0 3.0; 0.5 0.3]
        mpm = MatrixProjectionModel(A)
        @test mpm[1, 2] ≈ 3.0
        @test mpm[2, 1] ≈ 0.5
        @test Matrix(mpm) ≈ A
        @test length(mpm) == 4
    end

    @testset "Dimension validation" begin
        @test_throws DimensionMismatch MatrixProjectionModel(zeros(2, 3))
        @test_throws DimensionMismatch MatrixProjectionModel(zeros(2, 2), zeros(3, 3))
    end

    @testset "Trait types" begin
        @test LeslieMPM() isa AbstractMPMStructure
        @test LefkovitchMPM() isa AbstractMPMStructure
        @test DensityIndependent() isa AbstractDensityDependence
        @test DensityDependent() isa AbstractDensityDependence
        @test Deterministic() isa AbstractStochasticity
        @test StochasticKernelResampled() isa AbstractStochasticity
        @test StochasticParameterResampled() isa AbstractStochasticity
    end

    @testset "show methods" begin
        A = [0.0 3.0; 0.5 0.3]
        mpm = MatrixProjectionModel(A)
        buf = IOBuffer()
        show(buf, mpm)
        @test occursin("MatrixProjectionModel", String(take!(buf)))
    end
end
