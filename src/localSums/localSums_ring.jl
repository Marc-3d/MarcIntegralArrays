
"""
    In-place Local sums in a "ring ROI" (by sustracting a small rectangle from a big rectangle).
"""
function localsums!( 
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
        tmp = integralSum( 
            intArr, 
            Tuple(c) .- rad_in,
            Tuple(c) .+ rad_in,
            Tuple(c) .- rad_out, 
            Tuple(c) .+ rad_out,
            f
        )
        output[ c ] = op( output[c], tmp )
    end 
    return nothing
end