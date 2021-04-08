module FV3Kernels

export State

include("states.jl")

#=
export SillyType

struct SillyType{T}
    a :: T
    sin_a :: T
end

function SillyType(a)
    sin_a = sin(a)
    return SillyType(a, sin_a)
end
=#

end # module
