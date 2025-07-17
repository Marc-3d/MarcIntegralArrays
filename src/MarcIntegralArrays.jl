module MarcIntegralArrays

using ColorTypes

include("integralArrays/integralArrays.jl"); 
include("integralArrays/integralArraysL2.jl"); 
include("integralArrays/integralArraysL2M.jl"); 
include("integralArrays/IntegralArrays_extra.jl");

include("integralArrays/IntegralArrays_Colors.jl" )
include("integralArrays/IntegralVectorFields.jl")

# minmax isn't faster than max.( size(data), min.( 1, ... ) )
@inline  minmax( a::Int   , min::Int=1, max::Int=10 ) = a + (min-a)*Int(a<min) + (max-a)*Int(a>max)
@inline  minmax( a::Dims{2}, min::Int, max::Dims{2} ) = ( minmax( a[1], min, max[1] ), minmax( a[2], min, max[2] ) )
@inline  minmax( a::Dims{3}, min::Int, max::Dims{3} ) = ( minmax( a[1], min, max[1] ), minmax( a[2], min, max[2] ), minmax( a[3], min, max[3] ) )

@inline clipmin( a::Int, min::Int=1 ) = a + (min-a)*Int(a<min)
@inline clipmax( a::Int, max::Int=1 ) = a + (max-a)*Int(a>max) 

include("integralSums/IntegralSum.jl");
include("integralSums/IntegralSumN.jl"); 
include("integralSums/IntegralSum_ring.jl");
include("integralSums/IntegralSumN_ring.jl");

# local filters implemented with integral arrays
include("local_sum_operations.jl")
include("local_extrema.jl")
include("local_L2.jl");

include( "integralFilters/localSTD.jl" )

export integralArray, getIntegralArray


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
     Inspired by IntegralArrays.jl, I find it very handy to define a "getindex"
    that computes integral sums behind the scenes. 
    
    However, unlike IntegralArrays.jl, I am very much against the idea of
    using ClosedIntervals or any extenal data types for this. This is very
    common in Julia (specializing for highly specific types), which requires 
    developers to become familiar with the whole ecosystem of packages in a
    particular field.
    
    Instead, I opted for overriding the use of UnitRanges with my integral 
    arrays to express integral sums within the ROI defined by the rages. For
    instance, `intA[ 2:10, 4:8 ]` will compute the integral sum within the 
    rectangle from the top-left corner (2,4) to the bottom-right corner (10,8).
    
    In addition, I introduce the tuples as indices, to indicate the corners
    of the ROI of interest to compute sums. For instance, `intA[ (2,4), (10,8) ]`
    is equivalent to the above example, `intA[ 2:10, 4:8 ]`. 
=#

# intA[ 10 ] -> linear index 10
Base.@propagate_inbounds Base.getindex(A::IntegralArray, i::Int) = A.arr[i]

# intA[ 1, 3 ] -> cartesian index (1,3)
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,N}, i::Vararg{Int,N}) where {T,N} = A.arr[i...]

# intA[ 3, 10 ] -> integral sum within the ROI ( 3:10 ) 
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,1}, L::Int, R::Int ) where {T} = integralSum( A.arr, L, R ) 

# intA[ (4,5), (8,9) ] -> integral sum within the ROI ( 4:8, 5:9 ) 
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,N}, TLF::NTuple{N,Int}, BRB::NTuple{N,Int}) where {T,N} = integralSum( A.arr, TLF, BRB ) 

# intA[ 4:8 ] -> integral sum within the ROI ( 4:8 )
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,1}, ry::UnitRange{Int} ) where {T} = integralSum( A.arr, ry.start, ry.stop )

# intA[ 4:8, 4:9 ] -> integral sum within the ROI ( 4:8, 5:9 )
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,2}, ry::UnitRange{Int}, rx::UnitRange{Int}) where {T} = integralSum( A.arr, (ry.start,rx.start), (ry.stop,rx.stop) )

# intA[ 4:8, 4:9, 1:10 ] -> integral sum within the ROI ( 4:8, 5:9, 1:10 )
Base.@propagate_inbounds Base.getindex(A::IntegralArray{T,3}, ry::UnitRange{Int}, rx::UnitRange{Int}, rz::UnitRange{Int}) where {T} = integralSum( A.arr, (ry.start,rx.start,rz.start), (ry.stop,rx.stop,rz.stop ))

# multistep pipelines
# include("local_L2_segmentation.jl")

end # module MarcIntegralArrays

