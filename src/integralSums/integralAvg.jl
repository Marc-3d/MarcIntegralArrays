
"""
    Computes integral sum divided by the number of elements in the ROI: sum(I) / N 
"""
function integralAvg(
    intArr::AbstractArray{T,N}, 
    TL::Dims{N}, 
    BR::Dims{N},
    f::T=T(1)
) where {
    T,
    N
}
    return integralAvg_unsafe( 
        intArr, 
        minmax( TL, 1, size(intArr).-1 ), 
        minmax( BR, 1, size(intArr).-1 ),
        f
    )
end

function integralAvg(
    intArr::AbstractArray{T,1}, 
    L::Int, 
    R::Int,
    f::T=T(1)
)::T where {
    T
}
    return integralAvg_unsafe( 
        intArr, 
        minmax( L, 1, size(intArr,1)-1 ), 
        minmax( R, 1, size(intArr,1)-1 ),
        f
    )
end


""" 1D """
function integralAvg_unsafe( 
    intArr::AbstractArray{T,1}, 
    L::Int, 
    R::Int,
    f::T=T(1) 
)::T where {
    T<:AbstractFloat
}
    @inbounds begin  
        return f * ( intArr[ R+1 ] - intArr[ L ] ) / T( R - L + 1 ) 
    end
end

""" 1D, 2D, 3D """
function integralAvg_unsafe( 
    intArr::AbstractArray{T,N}, 
    TL::Dims{N}, 
    BR::Dims{N},
    f::T=T(1)
)::T where {
    T<:AbstractFloat, 
    N
}
    @inbounds begin  
        NN = prod( BR .- TL .+ 1 ) 
        return integralSum_unsafe( intArr, TL, BR, f ) / T( NN )
    end
end