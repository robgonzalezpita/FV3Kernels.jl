using Pkg

current_directory = @__DIR__

Pkg.activate(joinpath(current_directory, ".."))

using FV3Kernels, NCDatasets, Test


input_file_name = joinpath(current_directory, "../data/inputs/c_sw_12x24.nc")

input_file = NCDataset(input_file_name)

state = State(input_file)

@test state.isd == -2

# @show state