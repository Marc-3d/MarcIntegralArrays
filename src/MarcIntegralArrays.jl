module MarcIntegralArrays

include("IntegralArrays.jl"); 
include("IntegralArrays_extra.jl")
include("local_sum_operations.jl")
include("local_extrema.jl")
include("local_L2.jl");

export integralArray, getIntegralArray

greet() = print("Hello World!")

# Wrapping integral array in the "IntegralArray" type.
struct IntegralArray{T,N}
    arr::AbstractArray{T,N}
end

# Creating an integral array object from an input array (2D or 3D)
function getIntegralArray( input::AbstractArray{T,N}
                         ) where {T<:AbstractFloat,N}
    intA = IntegralArray{T,N}( zeros( T, size(input) .+ 1 ) );
    integralArray_unsafe!( intA.arr, input ); 
    return intA
end

# Taken from IntegralArrays.jl. These overwrites will probably be useful. 
Base.size(A::IntegralArray) = size(A.arr);
Base.axes(A::IntegralArray) = axes(A.arr);
Base.lastindex(A::IntegralArray) = A.arr[size(A)...]

#=
     Similar to IntegralArrays.jl, it is very handy to define a getindex
     that computed integral sums behind the scenes. unlike IntegralArrays.jl, 
     I am very much against the idea of using ClosedIntervals or any extenal
     data types or external notation.
=#

# intA[ 10 ] -> linear index 10
Base.@propagate_inbounds Base.getindex(A::IntegralArray, i::Int) = A.arr[i]

# intA[ 1, 3 ] -> cartesian index (1,3)
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,N}, i::Vararg{Int,N}) where {T,N} = A.arr[i...]

# intA[ (4,5), (8,9) ] -> integral sum within the ROI ( 4:8, 5:9 ) 
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,N}, TLF::NTuple{N,Int}, BRB::NTuple{N,Int}) where {T,N} = integralSum( A.arr, TLF, BRB ) 

# intA[ 4:8, 4:9 ] -> integral sum within the ROI ( 4:8, 5:9 )
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,2}, ry::UnitRange{Int}, rx::UnitRange{Int}) where {T} = integralSum( A.arr, (ry.start,rx.start), (ry.stop,rx.stop) )

# intA[ 4:8, 4:9, 1:10 ] -> integral sum within the ROI ( 4:8, 5:9, 1:10 )
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,3}, ry::UnitRange{Int}, rx::UnitRange{Int}, rz::UnitRange{Int}) where {T} = integralSum( A.arr, (ry.start,rx.start,rz.start), (ry.stop,rx.stop,rz.stop ))


end # module MarcIntegralArrays

