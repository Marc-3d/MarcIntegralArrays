""" 
    Computes in-place local sums * number of pixels in the rectangular local ROI.
"""
function localsumNs!( 
    output::AbstractArray{T,N},
    intArr::AbstractArray{T,N},
    rad::Dims{N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T,
    N
}
    @inbounds for c in CartesianIndices( output )
        tmp = integralSumN( 
            intArr, 
            Tuple(c) .- rad, 
            Tuple(c) .+ rad,
            f
        )
        output[ c ] = op( output[c], tmp ); 
    end 
    return nothing
end

function localsumNs!( 
    output::AbstractArray{T,N},
    intArr::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T,
    N
}
    @inbounds for c in CartesianIndices( output )
        tmp = integralSumN( 
            intArr, 
            Tuple(c) .- rad_in, 
            Tuple(c) .+ rad_in,
            Tuple(c) .- rad_out,
            Tuple(c) .+ rad_out,
            f
        )
        output[ c ] = op( output[c], tmp ); 
    end 
    return nothing
end
