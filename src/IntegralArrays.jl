
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
  accepts 1D, 2D or 3D arrays of 'eltype' Float32 or Float64. If your inputs 
  are integers, you should convert the inputs, e.g. by running:

     'intArr = integralArray( Float32.( input ) )'
"""
function integralArray( input::AbstractArray{T,N} 
                      ) where {T<:AbstractFloat,N}

    intArr = zeros( T, size(input) .+ 1 ); 
    integralArray_unsafe!( intArr, input); 
    return intArr
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
  next iteration, so we can cache this value. This is also true in 1D and 3D.
"""
function integralArray!( intArr::IntegralArray{T,N}, 
                         input::AbstractArray{T,N} 
                       ) where {T<:AbstractFloat,N}

     @assert size(intArr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr.arr, input ); 
     return nothing    
end

function integralArray!( intArr::AbstractArray{T,N}, 
                         input::AbstractArray{T,N} 
                       ) where {T<:AbstractFloat,N}

     @assert size(intArr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr, input ); 
     return nothing    
end
    
function integralArray_unsafe!( intArr::AbstractArray{T,1},
                                vec::AbstractArray{T,1} 
                              ) where {T<:AbstractFloat}
    cache = T(0)
    @inbounds for r in 2:size(vec,1)+1
          cache = cache + vec[r-1]
          intArr[r] = cache
    end
    return nothing
end

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

# The functions below accept a function-parameters, which is applied to the
# input data. this allows, for instance, to create "integral arrays squared",
# without having to create a squared input.

function integralArray!( intArr::IntegralArray{T,N}, 
                         input::AbstractArray{T,N},
                         f::Function=(x)::T->(x) 
                       ) where {T<:AbstractFloat,N}

     @assert size(intArr) == ( size(input) .+ 1 ); 
     integralArray_unsafe_f!( intArr.arr, input, f ); 
     return nothing    
end

function integralArray!( intArr::AbstractArray{T,N}, 
                         input::AbstractArray{T,N},
                         f::Function=(x)::T->(x) 
                       ) where {T<:AbstractFloat,N}

     @assert size(intArr) == ( size(input) .+ 1 ); 
     integralArray_unsafe_f!( intArr, input, f ); 
     return nothing    
end

function integralArray_f!( intArr::AbstractArray{T,N}, 
                           input::AbstractArray{T,N},
                           f::Function=(x)::T->(x)
                         ) where {T<:AbstractFloat,N}

     @assert size(intArr) == ( size(input) .+ 1 ); 
     integralArray_unsafe_f!( intArr, input, f ); 
     return nothing    
end

function integralArray_unsafe_f!( intArr::AbstractArray{T,1},
                                  vec::AbstractArray{T,1},
                                  f::Function=(x)::T->(x)
                                ) where {T<:AbstractFloat}
    cache = T(0)
    @inbounds for r in 2:size(vec,1)+1
          cache = cache + f(vec[r-1])
          intArr[r] = cache
    end
    return nothing
end

function integralArray_unsafe_f!( intArr::AbstractArray{T,2},
                                  img::AbstractArray{T,2},
                                  f::Function=(x)::T->(x)
                                ) where {T<:AbstractFloat}
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache = f(img[r-1,c-1]) + cache + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

function integralArray_unsafe_f!( intArr::AbstractArray{T,3},
                                  vol::AbstractArray{T,3},
                                  f::Function=(x)->(x)
                                ) where {T<:AbstractFloat}

    @inbounds for z in 2:size(vol,3)+1, c in 2:size(vol,2)+1
        cache = T(0)
        for r in 2:size(vol,1)+1
            cache += f(vol[r-1,c-1,z-1])
            intArr[r,c,z] = cache + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1];
        end
    end
    return nothing
end


#########


function test1()
   inp = rand( Float32, 100, 100 ); 
   TL = rand.( UnitRange.(1,size(inp)) ); 
   BR = rand.( UnitRange.(TL,size(inp)) ); 
   return sum( inp[TL[1]:BR[1],TL[2]:BR[2]] ) ==  MiA.integralSum( IA, TL, BR )
end
