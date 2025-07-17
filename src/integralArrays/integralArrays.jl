"""
    Base class
"""
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

"""
    Contructs an integral array from an input array of primitive numeric types (Ints or Floats). 
 
    IA = MarcIntegralArrays.IntegralArray( input )
 
    The last element in an integral array is the sum of all values in the input data. Thus, integral arrays are very sensitive to numer overflow. For example, a 50x50 (2500 pixels) image of type UInt8 will almost certainly have a total sum of intensities above 254 (the max value of UInt8). Therefore, integral arrays should always be "typed" with a higher bitdepth than the input data. To keep things simple, we use Float64 as a default type. You may try Float32 if you have good reason to belive 32 bits will not run into precission overflow.

    IA = MarcIntegralArrays.IntegralArray( input, Float32 ) 
"""
function IntegralArray( 
    inp::AbstractArray{R,N};
    T = Float64, 
    fun = (x)->(x)
) where {
    R<:Real,
    N
}
    IA = IntegralArray( T, size(inp) .+ 1 )
    integralArray!( IA, inp, fun )
    return IA
end

#=
   I like padding integral arrays with 0's at the beginning of each dim, because it helps for dealing with edge values when constructing the integral arrays. It also helps for computing integral sums adjacent to the top-left-front border.
  
  Because of this 1-element padding at the beginning of each dimension:
    intArr[2,2] = input[1,1]
=# 

"""
   IntegralArrays can be computed in O(N) by exploiting the fact that the cummulative sum at each position can be computed by reusing the cummulative sums previously computed in adjacent positions. For instance, each cummulative sum in a 2D  integral array can be computed by a combination of previously computed (left,top and top-left) cummulative sums in the integral array... plus the current element in the input data. More precisely:
    
  intArr[r,c] = input[r,c] + intArr[r-1,c] + intArr[r,c-1] - intArr[r-1,c-1]
		  
   Since we iterate row by row, intArr[r,c] becomes intArr[r-1,c] in the next iteration, so we can cache this value. This is also true in 1D and 3D.
"""
function integralArray!( 
    intArr::IntegralArray{T,N}, 
    input::AbstractArray{R,N}, 
    fun::Function=(v)->(v)
) where {
    T<:AbstractFloat,
    R<:Real,
    N
}
     @assert size(intArr.arr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr.arr, input, fun ); 
     return nothing    
end

function integralArray!( 
    intArr::AbstractArray{T,N}, 
    input::AbstractArray{R,N}, 
    fun::Function=(v)->(v)
) where {
    T<:AbstractFloat,
    R<:Real,
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
    vec::AbstractArray{R,1},
    fun::Function=(v)->(v)
) where {
    T<:AbstractFloat,
    R<:Real
}
    cache = T(0)
    @inbounds for r in 2:size(vec,1)+1
          cache += fun(T(vec[r-1]))
          intArr[r] = cache
    end
    return nothing
end

function integralArray_unsafe!( 
    intArr::AbstractArray{T,1},
    vec::AbstractArray{T,1},
    fun::Function=(v)->()
) where {
    T<:AbstractFloat
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
    img::AbstractArray{R,2},
    fun::Function=(v)->(v) 
) where {
    T<:AbstractFloat,
    R<:Real
}
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache += fun(T(img[r-1,c-1])) + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

function integralArray_unsafe!( 
    intArr::AbstractArray{T,2},
    img::AbstractArray{T,2},
    fun::Function=(v)->(v) 
) where {
    T<:AbstractFloat
}
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache += fun(img[r-1,c-1]) + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

#= 3D: It seems to be correct. I tested it with this code, and by changing parameters. 

    using MarcIntegralArrays
    x = rand( Float64, (10,10,10) )
    IA = marcIntegralArrays.IntegralArray( x ); 
    TLF = ( 1, 1, 1 ); 
    BRB = ( 3, 5, 2 );
    println( x[ UnitRange.( TLF, BRB )... ], "\t", IA[ TLF, BRB ] )
=#
function integralArray_unsafe!( 
    intArr::AbstractArray{T,3},
    vol::AbstractArray{R,3},
    fun::Function=(v)->(v) 
) where {
    T<:AbstractFloat,
    R<:Real
}
    @inbounds for z in 2:size(vol,3)+1, c in 2:size(vol,2)+1
        cache = T(0)
        for r in 2:size(vol,1)+1
            cache += fun(T(vol[r-1,c-1,z-1]))
            intArr[r,c,z] = cache + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1];
        end
    end
    return nothing
end

function integralArray_unsafe!( 
    intArr::AbstractArray{T,3},
    vol::AbstractArray{T,3},
    fun::Function=(v)->(v) 
) where {
    T<:AbstractFloat
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
