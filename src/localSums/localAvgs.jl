""" 
    input: 
        1.1) an array of numbers
        1.2) an array of grayscale values 

        2.1) a global "radius" 
        2.2) an array of potentially different "radii[ coord ]" for each element in the array

    output:
        1) an array{T} with the average of values around each pixel within a distance "radius" or "radii[ coord ]"
"""
function localAvgs( 
    input::AbstractArray{C,N}, 
    rad::Union{Dims{N},AbstractArray{Int,N}};
    T = Float64 
) where {
    C<:Union{Real,Color{<:Any,1}},
    N
}
    output = zeros( T, size(input));
    intA   = IntegralArray( input, T=T ); 
    localAvgs!( output, intA.arr, rad );
    return output;
end

""" 
    input: 
        1.1) an array of 3D-colors (RGB,HSL,YMC,etc)

        2.1) a global "radius" 
        2.2) an array of potentially different "radii[ coord ]" for each element in the array

    output:
        1) an array{T} with the average of values around each pixel within a distance "radius" or "radii[ coord ]"
"""
function localAvgs( 
    input::AbstractArray{C,N}, 
    rad::Union{Dims{N},AbstractArray{Int,N}};
    channels=collect(1:3),
    T = Float64 
) where {
    C<:Color{<:Any,3},
    N
}
    output = zeros( T, size(input)..., 3 );
    intA   = IntegralArray( input, T=T );
    out_v  = view( output, :, :, 1 )
    localAvgs!( out_v, intA.arr, rad );

    for i in 2:3
        c = channels[i]
        integralArray!( intA, input, channel=c ); 
        out_v = view( output, :, :, c )
        localAvgs!( out_v, intA.arr, rad )
    end

    colorType = C.name.wrapper
    return float2color3D( output, colorType );
end

# this is very slow
function float2color3D( inp::Array{T,3}, colorType=RGB ) where {T}
    oT  = colorType{T}
    out = zeros( oT, size(inp)[1:2] )
    @inbounds for c in CartesianIndices( out )
        c1, c2 = c.I
        out[c] = oT( inp[c1,c2,1], inp[c1,c2,2], inp[c1,c2,3]  )
    end;
    return out
end



""" 
    in-place local sums with a global radius
"""
function localAvgs!( 
    output::AbstractArray{T,N},
    intArr::AbstractArray{T,N},
    rad::Dims{N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(in)
) where {
    T,
    N
}
    @inbounds for c in CartesianIndices( output )
        tmp = integralAvg( 
            intArr, 
            Tuple(c) .- rad, 
            Tuple(c) .+ rad,
            f 
        )
        output[ c ] = op( output[c], tmp ); 
    end 
    return nothing
end

"""
    in-place local sums with a different radius for each element
"""
function localAvgs!( 
    output::AbstractArray{T,N},
    intArr::AbstractArray{T,N},
    rad::AbstractArray{Int,N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(in)
) where {
    T,
    N
}
    @inbounds for c in CartesianIndices( output )
        tmp = integralAvg( 
            intArr, 
            Tuple(c) .- rad[c], 
            Tuple(c) .+ rad[c],
            f 
        )
        output[ c ] = op( output[c], tmp ); 
    end 
    return nothing
end