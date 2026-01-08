function localAvgs( 
    input::AbstractArray{T,N},
    mask::AbstractArray{U,N},
    rad::Dims{N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T,
    U <: Real,
    N
}
    intAM = IntegralArrayM( input, mask )
    output = zeros( T, size(input) )
    localAvgs!( output, intAM, rad, f=f, op=op )
    return output
end

""" 
    Computes in-place local average ( local sums / number of non-masked pixels ) in the rectangular local ROI.
"""
function localAvgs!( 
    output::AbstractArray{T,N},
    intAM::IntegralArrayM{T,N},
    rad::Dims{N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T,
    N
}
    @inbounds for c in CartesianIndices( output )
        num = integralSum( 
            intAM.IA, 
            Tuple(c) .- rad, 
            Tuple(c) .+ rad,
            f 
        )
        den = integralSum( 
            intAM.IAM, 
            Tuple(c) .- rad, 
            Tuple(c) .+ rad,
            f 
        )
        tmp = num/den
        ( isnan( tmp ) || isinf( tmp ) ) && ( tmp = T(0); )
        output[ c ] = op( output[c], tmp ); 
    end 
    return nothing
end