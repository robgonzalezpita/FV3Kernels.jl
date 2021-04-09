#=
    Regression tests verify that the sw_driver program successfully ran for each datasize

=#

# include source file of the driver program
include("../src/c_sw.jl")

using Test, TOML, .program

configfile = TOML.parsefile("./test/data/inputs/inputs.toml")

# Iterate over the data sets in the config file
for (dataset, dataIODict) in configfile

    @testset "$dataset Regression " begin
        
        main(dataIODict)
        @test true

    end

end

return true