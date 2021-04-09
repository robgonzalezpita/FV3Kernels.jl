
# Unit Tests of functions within src/sw_core.jl

# using c_sw 

# using Test

function function_tests()
    @testset "c_sw function tests" begin
    # @test c_sw!()

    # @test divergence_corner()

    # @test d2a2c_vect()

    # @test edge_interpolate4()

    # @test fill2_4corners()

    # @test fill_4corners()

    end
end

@testset "All tests" begin
    function_tests()
end