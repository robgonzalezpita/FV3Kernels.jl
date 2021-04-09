#=

    Top level test suite, which runs unit tests, Pkg dependencies, regression & comparison tests.

=#

using Test

run(`cd ../`)

if VERSION >= v"1.6.0" 

    println("Starting All tests...")

    @testset "All tests" begin 
        include("test_c_sw.jl")
        
        # Both of the test suites below are a work in progress. 
        # include("unit_tests.jl")
        # include("PkgTests.jl")
    end

else

    println("Current version of Julia is not supported, please update to v1.6")

end
