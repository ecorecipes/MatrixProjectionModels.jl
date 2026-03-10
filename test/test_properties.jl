@testset "Properties" begin
    @testset "is_leslie" begin
        # Valid Leslie matrix
        A = [0.0 1.5 3.0; 0.8 0.0 0.0; 0.0 0.5 0.0]
        @test is_leslie(A) == true

        # Not Leslie: has diagonal element
        B = [0.0 1.5 3.0; 0.8 0.2 0.0; 0.0 0.5 0.0]
        @test is_leslie(B) == false

        # Not Leslie: missing sub-diagonal
        C = [0.0 1.5 3.0; 0.0 0.0 0.0; 0.0 0.5 0.0]
        @test is_leslie(C) == false

        # MatrixProjectionModel dispatch
        mpm = MatrixProjectionModel(A)
        @test is_leslie(mpm) == true
    end

    @testset "is_irreducible" begin
        # Irreducible: every stage reachable from every other
        A = [0.0 3.0; 0.5 0.0]
        @test is_irreducible(A) == true

        # Reducible: stage 2 can't reach stage 1
        B = [0.5 0.0; 0.3 0.5]
        @test is_irreducible(B) == false
    end

    @testset "is_primitive" begin
        # Primitive: irreducible + aperiodic
        A = [0.1 3.0; 0.5 0.1]
        @test is_primitive(A) == true

        # Imprimitive: irreducible but periodic (period 2)
        B = [0.0 3.0; 0.5 0.0]
        @test is_primitive(B) == false
    end

    @testset "is_ergodic" begin
        A = [0.1 3.0; 0.5 0.1]
        @test is_ergodic(A) == true
    end
end
