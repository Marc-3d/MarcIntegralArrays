"""
  Computes the sum of intensities within a rectangular ROI defined by its top-left[-front] (TL[F]) and bottom-right[-back] (BR[B]) corners. This is a 'safe' function (in contrast to `integralSum_unsafe`), which crops the input ROI to make sure it doesn't go out-of-bounds. 

By using the inlined function (minmax), the 'safe' code is on-par with the unsafe code in terms of performance. Therefore, you probably should stick to the 'safe' function, unless you are very desperate for performance.
"""
function integralSum( 
    intArr::IntegralArray{T,N}, 
    TL::Dims{N},
    BR::Dims{N},
    f::T=T(1)
) where {
    T,
    N
}
    return integralSum( intArr.arr, TL, BR, f )
end

function integralSum( 
    intArr::AbstractArray{T,N}, 
    TL::Dims{N},
    BR::Dims{N},
    f::T=T(1)
) where {
    T,
    N
}
    return integralSum_unsafe( 
      intArr, 
      max.( 1, min.( TL, size(intArr) .- 1  ) ),
      max.( 1, min.( BR, size(intArr) .- 1  ) ),
      f
    )
end

# 1D: using Int instead of of Dims{1}


function integralSum( 
  intArr::AbstractArray{T,1}, 
  L::Int,
  R::Int,
  f::T=T(1)
) where {
  T
}
  return integralSum_unsafe( 
    intArr, 
    max( 1, min( L, size(intArr,1) - 1 ) ), 
    max( 1, min( R, size(intArr,1) - 1 ) ),
    f
  )
end

############

function integralSum_unsafe( 
    intArr::IntegralArray{T,N}, 
    TL::Dims{N}, 
    BR::Dims{N},
    f::T=T(1)
) where {
    T<:AbstractFloat, 
    N
}
    return integralSum_unsafe( intArr.arr, TL, BR, f )
end


""" 1D integrals sums unsafe """
function integralSum_unsafe( 
    intArr::AbstractArray{T,1}, 
    L::Dims{1}, 
    R::Dims{1},
    f::T=T(1)
) where {
    T<:AbstractFloat
}
    return integralSum_unsafe( 
        intArr, 
        L[1], 
        R[1],
        f 
    )
end

function integralSum_unsafe( 
    intArr::AbstractArray{T,1}, 
    L::Int, 
    R::Int,
    f::T=T(1)
)::T where {
    T<:AbstractFloat
}
    @inbounds begin  
        return f * ( intArr[ R+1 ]  - intArr[ L ] )
    end
end

""" 2D integral sums unsafe """
function integralSum_unsafe( 
    intArr::AbstractArray{T,2}, 
    TL::Dims{2}, 
    BR::Dims{2},
    f::T=T(1)
) where {
    T<:AbstractFloat
}
    @inbounds begin  
        return f * ( intArr[TL[1],TL[2]] - intArr[BR[1]+1,TL[2]] - intArr[TL[1],BR[2]+1] + intArr[BR[1]+1,BR[2]+1] )
    end
end

""" 3D integral sums unsafe """
function integralSum_unsafe( 
    intArr::AbstractArray{T,3}, 
    TLF::Dims{3}, 
    BRB::Dims{3},
    f::T=T(1)
) where {
    T<:AbstractFloat
}
    @inbounds begin 
        return f * ( intArr[BRB[1]+1,BRB[2]+1,BRB[3]+1] - intArr[TLF[1],TLF[2],TLF[3]] - intArr[TLF[1],BRB[2]+1,BRB[3]+1] - intArr[BRB[1]+1,TLF[2],BRB[3]+1] - intArr[BRB[1]+1,BRB[2]+1,TLF[3]] + intArr[BRB[1]+1,TLF[2],TLF[3]] + intArr[TLF[1],BRB[2]+1,TLF[3]] + intArr[TLF[1],TLF[2],BRB[3]+1] )
    end    
end
