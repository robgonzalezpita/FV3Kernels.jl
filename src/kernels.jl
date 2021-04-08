using KernelAbstractions

@kernel function compute_vorticity!(vorticity, fx, fy, k)
    i, j = @index(Global, NTuple)
    vorticity[i, j, k] = fx[i, j-1] - fx[i, j] - fy[i-1, j] + fy[i, j]
end

#=
@inline δx⁻(i, j, f::AbstractArray) = @inbounds f[i, j] - f[i-1, j]
@inline δx⁺(i, j, f::AbstractArray) = @inbounds f[i+1, j] - f[i, j]

@inline δy(i, j, f::AbstractArray) = @inbounds f[i, j] - f[i, j-1]

@inline δx(i, j, f::Function, args...) = f(i, j, args...) - f(i-1, j, args...)
@inline δy(i, j, f::Function, args...) = f(i, j, args...) - f(i, j-1, args...)

@inline δ²x(i, j, f, args...) = δx(i, j, δx, f, args...)

fx(i, j, args...) = #
fy(i, j, args...) = #

@kernel function compute_vorticity!(vorticity, k, args...)
    i, j = @index(Global, NTuple)
    vorticity[i, j, k] = fx(i, j-1, args...) - fx(i, j, args...) - fy(i-1, j, args...) + fy(i, j, args...)
end
=#
