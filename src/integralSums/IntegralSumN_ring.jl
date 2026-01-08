
function integralSumN( 
  intArr::AbstractArray{T,N}, 
  TL1::Dims{N},
  BR1::Dims{N},
  TL2::Dims{N},
  BR2::Dims{N},
  f::T
) where {
  T,
  N
}
  return integralSumN( 
    intArr, 
    minmax.( TL1, 1, size(intArr).-1 ), 
    minmax.( BR1, 1, size(intArr).-1 ),
    minmax.( TL2, 1, size(intArr).-1 ), 
    minmax.( BR2, 1, size(intArr).-1 ),
    f
  )
end

function integralSumN( 
  intArr::AbstractArray{T,1}, 
  L1::Int,
  R1::Int,
  L2::Int,
  R2::Int,
  f::T
) where {
  T
}

  return integralSumN( 
    intArr, 
    minmax( L1, 1, size(intArr,1)-1 ), 
    minmax( R1, 1, size(intArr,1)-1 ),
    minmax( L2, 1, size(intArr,1)-1 ), 
    minmax( R2, 1, size(intArr,1)-1 ),
    f
  )
end


# A version that is aware of "N", the size of the ROI after bound-safe cropping,
# and combined the integralSum and N according to the "two-input" function "op".

""" 1D """
function integralSumN( 
  intArr::AbstractArray{T,1}, 
  L1::Int, 
  R1::Int, 
  L2::Int, 
  R2::Int,
  f::T
)::T where {
  T<:AbstractFloat
}
    @inbounds begin  
        return integralSum( intArr, L1, R1, L2, R2, f ) * T( R2 - L2 + 1 - R1 + L1 - 1 )
    end
end

#
function integralSumN( 
  intArr::AbstractArray{T,1}, 
  L1::Dims{1}, 
  R1::Dims{1}, 
  L2::Dims{1}, 
  R2::Dims{1},
  f::T
)::T where {
  T<:AbstractFloat
}
    return integralSumN( intArr, L1[1], R1[1], L2[1], R2[1], f )
end

""" 2D """
function integralSumN( 
  intArr::AbstractArray{T,2}, 
  TL1::Dims{2}, 
  BR1::Dims{2}, 
  TL2::Dims{2}, 
  BR2::Dims{2},
  f::T
)::T where {
  T<:AbstractFloat
}
    @inbounds begin  
        return integralSum( intArr, TL1, BR1, TL2, BR2, f ) *  T( prod( BR2 .- TL2 .+ 1 ) - prod( BR1 .- TL1 .+ 1) )
    end
end

""" 3D """
function integralSumN( 
  intArr::AbstractArray{T,3}, 
  TLF1::Dims{3}, 
  BRB1::Dims{3}, 
  TLF2::Dims{3}, 
  BRB2::Dims{3},
  f::T
)::T where {
  T<:AbstractFloat
}
    @inbounds begin 
        return integralSum( intArr, TLF1, BRB1, TLF2, BRB2, f ) * T( prod( BRB2 .- TLF2 .+ 1 ) - prod( BRB1 .- TLF1 .+ 1 ) )
    end    
end

