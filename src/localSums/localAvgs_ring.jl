
function localAvgs( 
    input::AbstractArray{R,N},
    rad1::Dims{N},
    rad2::Dims{N};
    T=Float64,
    f=T(1),
    op::Function=(out,in)->(out+in)
) where {
    R<:Real,
    N
}
    intA = IntegralArray( input, T )
    output = zeros( T, size(input) )
    localAvgs!( output, intA.arr, rad1, rad2, f=f, op=op )
    return output
end


function localAvgs!( 
    output::AbstractArray{T,N},
    intArr::AbstractArray{T,N},
    rad1::Dims{N},
    rad2::Dims{N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T,
    N
}
    @inbounds for c in CartesianIndices( output )
        tmp = integralAvg( 
            intArr, 
            Tuple(c) .- rad1, 
            Tuple(c) .+ rad1,
            Tuple(c) .- rad2, 
            Tuple(c) .+ rad2,
            f 
        )
        output[ c ] = op( output[c], tmp ); 
    end 
    return nothing
end
