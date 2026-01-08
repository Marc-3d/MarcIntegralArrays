"""
    Integral arrays for integral masked sums and averages
"""
mutable struct IntegralArrayM{T,N}
    IA::IntegralArray{T,N}
    IAM::IntegralArray{T,N}
end

# default constructor
function IntegralArrayM( 
    type, 
    size::Dims{N} 
) where {
    N
}
    return IntegralArrayM{type,N}( 
        IntegralArray( type, size ),
        IntegralArray( type, size ),
    )
end

# contructor from input data
function IntegralArrayM( 
    inp::AbstractArray{I,N},
    mask::AbstractArray{U,N},
    T=Float64
) where {
    I<:Union{Real,Color{<:Any,1}},
    U<:Real,
    N
}
    IAM = IntegralArrayM( T, size(inp) .+ 1 )
    integralArrayM!( IAM, inp, mask )

    return IAM
end

# in-place computation of the two integral arrays 
function integralArrayM!( 
    IAM::IntegralArrayM{T,N}, 
    inp::AbstractArray{I,N},
    mask::AbstractArray{U,N}
) where {
    T<:AbstractFloat,
    I<:Union{Real,Color{<:Any,1}},
    U<:Real,
    N
}
    # TODO: integral array that multiplies the mask within the intgral array loop?
    integralArray!( IAM.IA , inp .* mask )
    integralArray!( IAM.IAM, mask )
    return nothing
end
