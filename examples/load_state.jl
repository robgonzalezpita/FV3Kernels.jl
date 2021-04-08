using FV3Kernels

input_file_name = "something.nc"

state = State(input_file_name)

# @test state.isd == 7
