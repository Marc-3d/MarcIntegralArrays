
"""
    This class is meant two store two integral arrays for the desired input array. Namely, it hold a "standard integral array", and a "square integral array". Both of these are required to compute local L2s.
"""
mutable struct IntegralArraysL2{T,N}
    IA::IntegralArray{T,N}
    IA2::IntegralArray{T,N}
end

# default constructor
function IntegralArraysL2( 
    T::Type, 
    size::Dims{N} 
) where {
    N
}
    return IntegralArraysL2{T,N}( 
        IntegralArray( T, size ), 
        IntegralArray( T, size ) 
    )
end

# contructor from input data
function IntegralArraysL2( 
    inp::AbstractArray{I,N},
    T=Float64
) where {
    I<:Real,
    N
}
    # initializing the integral arrays within IAL2
    IAL2 = IntegralArraysL2( T, size(inp) .+ 1 )
    # populating the integral arrays within IAL2
    integralArraysL2!( IAL2, inp )

    return IAL2
end

# in-place computation of the two integral arrays 
function integralArraysL2!( 
    IAL2::IntegralArraysL2{T,N}, 
    inp::AbstractArray{I,N} 
) where {
    T<:AbstractFloat,
    I<:Real,
    N
}
    integralArray!( IAL2.IA , inp )
    integralArray!( IAL2.IA2, inp, (x)->(T(x)^2) )
    return nothing
end
