""" 
    input: 
        1) arrays of numbers and grayscale values 
        2) a global "radius" or an array of pontentially different "radii[ coord ]" for each element in the array
    output:
        1) an array with the sum of values around each pixel within a distance "radius" or "radii[ coord ]"
"""
function localsums( 
    input::AbstractArray{C,N}, 
    rad::Union{Dims{N},AbstractArray{Int,N}};
    T = Float64 
) where {
    C<:Union{Real,Color{<:Any,1}},
    N
}
    output = zeros( T, size(input));
    intA   = IntegralArray( input, T=T ); 
    localsums!( output, intA.arr, rad );
    return output;
end

""" arrays of 3D colors (RGB, YMC, HSL) """
function localsums( 
    input::AbstractArray{C,N}, 
    rad::Dims{N};
    channels=collect(1:3),
    T = Float64 
) where {
    C<:Color{<:Any,3},
    N
}
    output = zeros( T, size(input)..., 3 );
    intA   = IntegralArray( input, T=T );
    out_v  = view( output, :, :, 1 )
    localsums!( out_v, intA.arr, rad );

    for i in 2:3
        c = channels[i]
        integralArray!( intA, input, channel=c ); 
        out_v = view( output, :, :, c )
        localsums!( out_v, intA.arr, rad )
    end

    # this is averaging... should go into "localAvgs"
    #NN = prod( 2 .* rad .+ 1 ) 
    #output ./= NN

    return float2RGB( output );
end

function localsums( 
    input::AbstractArray{C,N}, 
    rad::AbstractArray{Int,N};
    channels=collect(1:3),
    T = Float64 
) where {
    C<:Color{<:Any,3},
    N
}
    output = zeros( T, size(input)..., 3 );
    intA   = IntegralArray( input, T=T );
    out_v  = view( output, :, :, 1 )
    localsums!( out_v, intA.arr, rad );

    for i in 2:3
        c = channels[i]
        integralArray!( intA, input, channel=c ); 
        out_v  = view( output, :, :, c )
        localsums!( out_v, intA.arr, rad )
    end

    # TOOD: hard-corded for 2D at the moment
    # this is averaging... should go into "localAvgs"
    #for c in CartesianIndices( input )
    #    output[Tuple(c)...,:] ./= prod( 2 .* (rad[c],rad[c]) .+ 1 )
    #end

    return float2RGB( output );
end

function float2RGB( inp::Array{T,3} ) where {T}
    out = zeros( RGB{T}, size(inp)[1:2] )
    for c in CartesianIndices( out )
       out[c] = RGB{T}( inp[Tuple(c)...,:]... )
    end;
    return out
end



""" 
    in-place local sums with a global radius
"""
function localsums!( 
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
        tmp = integralSum( 
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
function localsums!( 
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
        tmp = integralSum( 
            intArr, 
            Tuple(c) .- rad[c], 
            Tuple(c) .+ rad[c],
            f 
        )
        output[ c ] = op( output[c], tmp ); 
    end 
    return nothing
end