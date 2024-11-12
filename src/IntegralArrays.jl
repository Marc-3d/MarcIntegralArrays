
#=
   INTEGRAL ARRAYS : CREATING THEM AND COMPUTING SUMS WITH THEM

   I like padding integral arrays with 0's at the beginning of each dim, 
  because it helps for dealing with edge values when constructing the  
  integral arrays. It also helps for computing integral sums adjacent to
  the top-left-front border.
  
  Because of this 1-element padding at the beginning of each dimension:
    intArr[2,2] = input[1,1]
=# 

"""
   This function returns the integral array of the input array. The function 
  accepts 2D or 3D arrays of 'eltype' Float32 or Float64. If your inputs are
  integers, you should convert the inputs, e.g. by running:

     'intArr = integralArray( Float32.( input ) )'
"""
function integralArray( input::AbstractArray{T,N} 
                      ) where {T<:AbstractFloat,N}

    intArr = zeros( T, size(input) .+ 1 ); 
    integralArray_unsafe!( intArr, input); 
    return intArr
end

function integralArray!( intArr::AbstractArray{T,N}, 
                         input::AbstractArray{T,N} 
                       ) where {T<:AbstractFloat,N}

     @assert size(intArr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr, input ); 
     return nothing    
end
    
"""
   IntegralArrays can be computed in O(N) by exploiting the fact that the
  cummulative sum at each position can be computed by reusing the cummulative
  sums in previous adjacent positions. For instance, each cummulative sum 
  in a 2D  integral array can be computed by a combination of previously 
  computed sums in the integral array (plus one element in the input data). 

  To be more precise:
    
  intArr[r,c] = input[r,c] + intArr[r-1,c] + intArr[r,c-1] - intArr[r-1,c-1]
		  
   Since we iterate row by row, intArr[r,c] becomes intArr[r-1,c] in the
  next iteration, so we can cache this value. This is also true in 3D.
"""
function integralArray_unsafe!( intArr::AbstractArray{T,2},
                                img::AbstractArray{T,2} 
                              ) where {T<:AbstractFloat}
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache = img[r-1,c-1] + cache + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

function integralArray_unsafe!( intArr::AbstractArray{T,3},
                                vol::AbstractArray{T,3}
                              ) where {T<:AbstractFloat}

    @inbounds for z in 2:size(vol,3)+1, c in 2:size(vol,2)+1
        cache = T(0)
        for r in 2:size(vol,1)+1
            cache += vol[r-1,c-1,z-1]
            intArr[r,c,z] = cache + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1];
        end
    end
    return nothing
end
    
# INTEGRAL SUMS
# NOTE: All benchmarks below are computed on a machine with an Intel Core i7-10875H.

@inline minmax( a, min=1, max=10 ) = a + (min-a)*(a<min) + (a>max)*(max-a)
@inline clipmin( a, min=1 ) = a + (min-a)*(a<min)
@inline clipmax( a, max=1 ) = a + (a>max)*(max-a)

"""
   Computes the sum of intensities within a rectangular ROI defined by its 
  top-left[-front] (TL[F]) and bottom-right[-back] (BR[B]) corners. The 'safe'
  function, 'integralSum', ensures that the input corners (TLF,BRB) are within
  bounds, which makes its slower than the unsafe version, 'integralSum_unsafe'.
  You probably should stick to the 'safe' function, unless you are desperate for
  performance.
"""

# 🚀 @btime integralSum( $intArr, $TL, $BR ) # in 2D
#      5.500 ns (0 allocations: 0 bytes)
function integralSum( intArr::AbstractArray{T,N}, 
                      TL::Dims{N},
                      BR::Dims{N} 
                    ) where {T,N}

    return integralSum_unsafe( intArr, 
                               minmax.( TL, 1, size(intArr).-1 ), 
                               minmax.( BR, 1, size(intArr).-1 ) );
end

# 🚀 @btime integralSum_unsafe( $intArr, $TL, $BR )
#      2.400 ns (0 allocations: 0 bytes)
function integralSum_unsafe( intArr::AbstractArray{T,2}, 
                             TL::Dims{2}, 
                             BR::Dims{2} 
                           ) where {T<:AbstractFloat}
    @inbounds begin  
        return intArr[TL[1], TL[2] ] - intArr[BR[1]+1, TL[2] ] - intArr[TL[1],BR[2]+1] + intArr[BR[1]+1,BR[2]+1] 
    end
end
 
# 🚀 @btime integralSum_unsafe( $intA, $TLF, $BRB )
#      4.900 ns (0 allocations: 0 bytes)
function integralSum_unsafe( intArr::AbstractArray{T,3}, 
                             TLF::Dims{3}, 
                             BRB::Dims{3} 
                           ) where {T<:AbstractFloat}
    @inbounds begin 
        return intArr[BRB[1]+1,BRB[2]+1,BRB[3]+1] - intArr[TLF[1],TLF[2],TLF[3]] - intArr[TLF[1],BRB[2]+1,BRB[3]+1] - intArr[BRB[1]+1,TLF[2],BRB[3]+1] - intArr[BRB[1]+1,BRB[2]+1,TLF[3]] + intArr[BRB[1]+1,TLF[2],TLF[3]] + intArr[TLF[1],BRB[2]+1,TLF[3]] + intArr[TLF[1],TLF[2],BRB[3]+1]
    end    
end

"""
  Specialized functions to compute averages (sums divided by number of elements).
"""
function integralAverage( intArr::AbstractArray{T,N}, TL::Dims{N}, BR::Dims{N} ) where {T,N}
    return integralSumN_unsafe( intArr, minmax.( TL, 1, size(intArr).-1 ), minmax.( BR, 1, size(intArr).-1 ), op=/ );
end

#= 
  The functions below accept another function as the last argument, 
  'op', which takes as inputs the integralSum and the size of the 
  ROI. By setting 'op=/', the function will return the integral sum
  divided by the number of element in ROI, which is the average value
  in ROI. 'op=*' is used 'local_L2'. 
=#

function integralSumN( intArr::AbstractArray{T,N},
                       TL::Dims{N}, 
                       BR::Dims{N}, 
                       op::Function=/ 
                     )::T where {T,N}

    return integralSumN_unsafe( intArr, 
                                minmax.( TL, 1, size(intArr).-1 ), 
                                minmax.( BR, 1, size(intArr).-1 ), 
                                op );
end

# 🚀 @btime integralSumN_unsafe( $intA, $TL, $BR, $op )
#     2.900 ns (0 allocations: 0 bytes)
function integralSumN_unsafe( intArr::AbstractArray{T,2}, 
                              TL::Dims{2}, 
                              BR::Dims{2}, 
                              op::Function=/ 
                            )::T where {T<:AbstractFloat}
    @inbounds begin  
        return op( ( intArr[TL[1], TL[2] ] - intArr[BR[1]+1, TL[2] ] - intArr[TL[1],BR[2]+1] + intArr[BR[1]+1,BR[2]+1] ), T( prod( BR .- TL .+1 ) ) )
    end
end

# 🚀 @btime integralAverage_unsafe( $intA, $TLF, $BRB, $op )
#      5.700 ns (0 allocations: 0 bytes)
function integralSumN_unsafe( intArr::AbstractArray{T,3}, 
                              TLF::Dims{3}, 
                              BRB::Dims{3}, 
                              op::Function=/ 
                            ) where {T<:Real}
    @inbounds begin 
        return op( ( intArr[BRB[1]+1,BRB[2]+1,BRB[3]+1] - intArr[TLF[1],TLF[2],TLF[3]] - intArr[TLF[1],BRB[2]+1,BRB[3]+1] - intArr[BRB[1]+1,TLF[2],BRB[3]+1] - intArr[BRB[1]+1,BRB[2]+1,TLF[3]] + intArr[BRB[1]+1,TLF[2],TLF[3]] + intArr[TLF[1],BRB[2]+1,TLF[3]] + intArr[TLF[1],TLF[2],BRB[3]+1] ), prod( BRB .- TLF .+ 1 ) )
    end    
end





function test1()
   inp = rand( Float32, 100, 100 ); 
   TL = rand.( UnitRange.(1,size(inp)) ); 
   BR = rand.( UnitRange.(TL,size(inp)) ); 
   return sum( inp[TL[1]:BR[1],TL[2]:BR[2]] ) ==  MiA.integralSum( IA, TL, BR )
end
