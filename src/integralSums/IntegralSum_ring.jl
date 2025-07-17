"""
    This function accepts two radii, and computes the sum in a ring by substracting the sums in small ROI to the sums in a bigger ROI. The first set of ROI coordinates (TL1 and BR1) define the smaller ROI, and the TL2 and BR2 define the larger ROI. It basically involves computing two integral sums, and it is about twice as slow as computing a single integral sum.
"""
function integralSum( 
  intArr::AbstractArray{T,N}, 
  TL1::Dims{N},
  BR1::Dims{N},
  TL2::Dims{N},
  BR2::Dims{N},
  f::T=T(1)
) where {
  T,
  N
}
    return integralSum( 
        intArr, 
        TL2, 
        BR2,
        f
    ) - integralSum(
        intArr, 
        TL1, 
        BR1,
        f  
    )
end

# 1D: using Int instead of of Dims{1}
function integralSum(
  intArr::AbstractArray{T,1}, 
  L1::Int,
  R1::Int,
  L2::Int,
  R2::Int,
  f::T=T(1) 
) where {
  T
}
    return integralSum( 
        intArr, 
        (L2,), 
        (R2,),
        f
    ) - integralSum(
        intArr,
        (L1,), 
        (R1,),
        f
    )
end
