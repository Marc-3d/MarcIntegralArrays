
struct IntegralArray{T,N}
    arr::AbstractArray{T,N}
end

# default constructor
function IntegralArray( 
    T::Type, 
    size::Dims{N} 
) where {
    N
}
    return IntegralArray{T,N}( 
        zeros( T, size )
    )
end

# contructor from input data
function IntegralArray( 
    inp::AbstractArray{I,N},
    T = Float64
) where {
    I<:Real,
    N
}
    println(".....")
    IA = IntegralArray( T, size(inp) .+ 1 )
    integralArray!( IA, inp, (x)->(T(x)) )

    return IA
end

#=
   INTEGRAL ARRAYS : CREATING THEM AND COMPUTING SUMS WITH THEM

   I like padding integral arrays with 0's at the beginning of each dim, because it helps for dealing with edge values when constructing the integral arrays. It also helps for computing integral sums adjacent to the top-left-front border.
  
  Because of this 1-element padding at the beginning of each dimension:
    intArr[2,2] = input[1,1]
=# 

"""
   This function returns the integral array of the input array. The function accepts 1D, 2D or 3D input arrays of 'eltype' Float32 or Float64. If your inputs are integers, you should convert the inputs, e.g. by running:

     'intArr = integralArray( Float32.( input ) )'
"""
function integralArray( 
    input::AbstractArray{I,N},
    T::Type=Float64
) where {
    I<:Real,
    N
}
    intArr = zeros( T, size(input) .+ 1 ); 
    integralArray_unsafe!( intArr, input ); 
    return intArr
end

"""
   IntegralArrays can be computed in O(N) by exploiting the fact that the cummulative sum at each position can be computed by reusing the cummulative sums previously computed in adjacent positions. For instance, each cummulative sum in a 2D  integral array can be computed by a combination of previously computed (left,top and top-left) cummulative sums in the integral array... plus the current element in the input data. More precisely:
    
  intArr[r,c] = input[r,c] + intArr[r-1,c] + intArr[r,c-1] - intArr[r-1,c-1]
		  
   Since we iterate row by row, intArr[r,c] becomes intArr[r-1,c] in the next iteration, so we can cache this value. This is also true in 1D and 3D.
"""
function integralArray!( 
    intArr::IntegralArray{T,N}, 
    input::AbstractArray{I,N}, 
    fun::Function=(v)->(T(v))
) where {
    T<:AbstractFloat,
    I<:Real,
    N
}
     @assert size(intArr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr.arr, input, fun ); 
     return nothing    
end

function integralArray!( 
    intArr::AbstractArray{T,N}, 
    input::AbstractArray{I,N}, 
    fun::Function=(v)->(T(v))
) where {
    T<:AbstractFloat,
    I<:Real,
    N
}
     @assert size(intArr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr, input, fun ); 
     return nothing    
end

#############################

#= All tests bellow share the same total number of array elements: 10000.

using BenchmarkTools, MarcIntegralArrays

Ny = 10000; T = Float32; inp = zeros( Int32, Ny ); IA = zeros( T, Ny+1 );
@btime MarcIntegralArrays.integralArray_unsafe!( $IA, $inp )
  ~ 13 μs (AMD EPYC 7453)

Ny, Nx = 100, 100; T = Float32; inp = zeros( Int32, Ny, Nx ); IA = zeros( T, Ny+1, Nx+1 );
@btime MarcIntegralArrays.integralArray_unsafe!( $IA, $inp )
  ~ 15 μs (AMD EPYC 7453)

Ny, Nx, Nz = 10, 10, 100; T = Float32; inp = zeros( Int32, Ny, Nx, Nz ); IA = zeros( T, Ny+1, Nx+1, Nz+1 );
@btime MarcIntegralArrays.integralArray_unsafe!( $IA, $inp )
  ~ 12 μs (AMD EPYC 7453)
=#

# 1D
function integralArray_unsafe!( 
    intArr::AbstractArray{T,1},
    vec::AbstractArray{I,1},
    fun::Function=(v)->(T(v))
) where {
    T<:AbstractFloat,
    I<:Real
}
    cache = T(0)
    @inbounds for r in 2:size(vec,1)+1
          cache += fun(vec[r-1])
          intArr[r] = cache
    end
    return nothing
end

# 2D
function integralArray_unsafe!( 
    intArr::AbstractArray{T,2},
    img::AbstractArray{I,2},
    fun::Function=(v)->(T(v)) 
) where {
    T<:AbstractFloat,
    I<:Real
}
    #println( typeof(10), typeof( fun(10)) )
    #println( T )
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache += fun(img[r-1,c-1]) + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

# 3D: might be wrong
function integralArray_unsafe!( 
    intArr::AbstractArray{T,3},
    vol::AbstractArray{I,3},
    fun::Function=(v)->(T(v)) 
) where {
    T<:AbstractFloat,
    I<:Real
}
    @inbounds for z in 2:size(vol,3)+1, c in 2:size(vol,2)+1
        cache = T(0)
        for r in 2:size(vol,1)+1
            cache += fun(vol[r-1,c-1,z-1])
            intArr[r,c,z] = cache + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1];
        end
    end
    return nothing
end