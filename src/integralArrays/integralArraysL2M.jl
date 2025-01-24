
"""
    Integral arrays for integral masked squarred errors
"""
mutable struct IntegralArraysL2M{T,N}
    IA::IntegralArray{T,N}
    IA2::IntegralArray{T,N}
    IAM::IntegralArray{T,N}
end

# default constructor
function IntegralArraysL2M( 
    type, 
    size::Dims{N} 
) where {
    N
}
    return IntegralArraysL2M{type,N}( 
        IntegralArray( type, size ), 
        IntegralArray( type, size ),
        IntegralArray( type, size ),
    )
end

# contructor from input data
function IntegralArraysL2M( 
    inp::AbstractArray{I,N},
    mask::AbstractArray{U,N},
    T=Float64
) where {
    I<:Real,
    U<:Integer,
    N
}
    # initializing the integral arrays within IAL2
    IAL2 = IntegralArraysL2M( T, size(inp) .+ 1 )
    # populating the integral arrays within IAL2
    integralArraysL2M!( IAL2, inp, mask )

    return IAL2
end

# in-place computation of the two integral arrays 
function integralArraysL2M!( 
    IAL2M::IntegralArraysL2M{T,N}, 
    inp::AbstractArray{I,N},
    mask::AbstractArray{U,N}
) where {
    T<:AbstractFloat,
    I<:Real,
    U<:Integer,
    N
}
    integralArray!( IAL2M.IA , inp )
    integralArray!( IAL2M.IA2, inp, (x)->(T(x)^2) )
    integralArray!( IAL2M.IAM, mask )
    return nothing
end
