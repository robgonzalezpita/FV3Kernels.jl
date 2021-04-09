#=
      Driver program to run the c_sw kernel

        This program acheives the goals of:
            * reading a NetCDF file
            * printing statistics on the input & output state in ./test/test_output
            * running the c_sw kernel located in sw_core.jl
            * creating a new NetCDF file with the data after the kernel runs

=#

module program

include("./sw_core.jl")

using TOML, Printf
using NCDatasets, OffsetArrays
using .sw_core_mod


export main


    function main(input::Dict)

        # define vars
        datasize = input["name"]
        nc_input_file = input["input_file"]
        print_output_file = input["test_out_file"]
        nc_output_file = input["output_file"]


    # Define the dataset using NCDatasets pkg from each datafile 
    ds = NCDataset(nc_input_file, "a")

    # Open the export file with write/create/truncate privileges ("w")
    io = open(print_output_file,"w")

    # Assign Variables from the NetCDF file to a Julia Struct
    #
    # Could be written 
    #
    # current_state = State(ds)
    # different_state = State(different_ds)
    assign_variables(ds)


    # Close the NetCDF dataset
    close(ds)


    # Print input state 
    print_state("Input State", current_state, io)


    # Call the Julia kernel
    println("Run the kernel $datasize : ")
    @time for k = 1 : current_state.npz
        c_sw!(current_state, k)
    end


    # Print output state 
    print_state("Output State", current_state, io)


    # Write a new NetCDF file
    println("Write NetCDF file $datasize : ") 
    @time write_state(nc_output_file, current_state)


    # Close the IO file
    close(io)

    return true


    end # function main


end # module program
