"""
     For each position in the input data, computes the local sum 
    of values within a rectangular  ROI of half-size "rad". This 
    is a 'safe' function, meaning that it corrects for out-of-
    bounds coordinates. Thus, this function can compute local
    sums at positions that are less than "rad-distance" away from
    the input borders. However, the ROI's at these positions are 
    cropped to remove out-of-bounds regions, which means that the
    sums around the edges are artifically low. 
"""

# inp = rand( Float32, 100, 100 ); r = (10,10)
# @btime out = localsums( $inp, $r )
#    74.600 μs (4 allocations: 79.03 KiB); 
function localsums( input::AbstractArray{T,N}, 
                    rad::Dims{N} ) where {T<:AbstractFloat,N}

    output = zeros( T, size(input));
    intA   = integralArray( input ); 
    localsums!( intA, output, rad );
    return output;
end

function localsums!( intArr::AbstractArray{T,N},
                     output::AbstractArray{T,N},
	             rad::Dims{N} ) where {T,N}

    @inbounds for c in CartesianIndices( output )
        output[ c ] = integralSum( intArr, Tuple(c) .- rad, Tuple(c) .+ rad )
    end 
    return nothing
end

"""
     Computes the local sums within a rectangular ROI of half-size
    'rad' around each position of the input array. This function
    doesn't perform any bounds checks or correct for out-of-bounds
    coordinates, so it cannot safely deal with positions closer than 
    "rad-distance" to the edges. Thus, the function simply ignores
    these positions. Arguably, ignoring them can be better than computing
    artifiically lower sums (see 'localsums'). On the other hand, this
    function is faster than the 'safe' alternative.
"""

# r = (10,10); inp = rand( Float32, ( 100, 100 ) .+ 2 .* r ); 
#  @btime out = localsums_unsafe( $inp, $r )
#    44.300 μs (4 allocations: 113.59 KiB)
function localsums_unsafe( input::AbstractArray{T,N}, 
                           rad::Dims{N} ) where {T<:AbstractFloat,N}

    output = zeros( T, size(input));
    intA   = integralArray( input ); 
    localsums_unsafe!( intA, output, rad );
    return output;
end

# 2D
function localsums_unsafe!( intArr::AbstractArray{T,2}, 
                            output::AbstractArray{T,2}, 
                            rad::Dims{2} ) where {T}

    out_sz = size( output )
    out_0  = (1,1)  .+ rad
    out_1  = out_sz .- rad
    ROI_   = UnitRange.( out_0, out_1 ); 

    @inbounds for c in CartesianIndices( ROI_ )
        output[ c ] = integralSum_unsafe( intArr, Tuple(c) .- rad, Tuple(c) .+ rad )
    end
    return nothing
end

# 3D
function localsums_unsafe!( intArr::AbstractArray{T,3}, 
                            output::AbstractArray{T,3}, 
                            rad::Dims{3} ) where {T}

    out_sz = size( output )
    out_0  = (1,1,1) .+ rad
    out_1  = out_sz  .- rad
    ROI_   = UnitRange.( out_0, out_1 ); 

    @inbounds for c in CartesianIndices( ROI_ )
        output[ c ] = integralSums_unsafe( intArr, Tuple(c) .- rad, Tuple(c) .+ rad )
    end
    return nothing
end

"""
     Computes the local average around each position based on
    a rectangular ROI with half-size "rad". This function 
    uses 'safe' integralSums, so it will compute the average
    for position close the edges of the input. Unlike 'localsums',
    local averages are not biased when out-of-bounds ROIs are
    cropped, since the average value is relative to the number
    of elements in the ROI, not it's size.
"""
function localavg( input::AbstractArray{T,N}, 
                   rad::Dims{N} ) where {T<:AbstractFloat,N}

    output = zeros( T, size(input));
    intA   = integralArray( input ); 
    localavg!( intA, output, rad );
    return output;
end

function localavg!( intArr::AbstractArray{T,2},
                    output::AbstractArray{T,2},
	            rad::Dims{2} ) where {T}

    @inbounds for c in CartesianIndices( output )
        output[ c ] = integralAverage( intArr, Tuple(c) .- rad, Tuple(c) .+ rad )
    end 
    return nothing
end

"""
    This function is the holy grail for in-place localsums, as 
   it exposes parameters that control how the local sums are
   combined with the output array. In other words, it allows
   to do in-place addition, substraction, multiplication, etc of
   local sums. This is clearly useful in 'local_L2'.

    If you simply want to set each pixel to the local sum, either
   use 'localsums_unsafe' or pass the following op:
   
      localsums_unsfae_op!( intA, output, rad, (a,b)->(b) )

   It doesn't seem to add any significant overhead :)
"""

# r = (10,10); inp = rand( Float32, 100 + 2*r[1], 100 + 2*r[2] ); 
# @btime out = localsums_unsafe( $inp, $r )  
#   47.300 μs (4 allocations: 113.59 KiB)
function localsums_unsafe_f_op( input::AbstractArray{T,N}, 
                                rad::Dims{N},
                                f::T=T(1),
                                op::Function=+
                              ) where {T<:AbstractFloat,N}
    output = zeros( T, size(input));
    intA   = integralArray( input ); 
    localsums_unsafe_f_op!( intA, output, rad, f, op );
    return output;
end

# 2D
function localsums_unsafe_f_op!( intArr::AbstractArray{T,2}, 
                                 output::AbstractArray{T,2}, 
                                 rad::Dims{2}, 
                                 f::T=T(1),
                                 op::Function=+ ) where {T<:AbstractFloat}
    out_sz = size( output )
    out_0  = (1,1)  .+ rad
    out_1  = out_sz .- rad
    ROI_   = UnitRange.( out_0, out_1 ); 
    WH     = 2 .* rad

    @inbounds for c in CartesianIndices( ROI_ )
        output[ c ] = op( output[c], f * integralSum_unsafe( intArr, Tuple(c) .- rad, Tuple(c) .+ rad ) )
    end
    return nothing
end

# 3D
function localsums_unsafe_f_op!( intArr::AbstractArray{T,3}, 
                                 out::AbstractArray{T,3}, 
                                 rad::Dims{3},
                                 f::T=T(1),
                                 op::Function=+ )  where {T<:AbstractFloat}

    out_sz = size( output )
    out_0  = (1,1,1) .+ rad
    out_1  = out_sz  .- rad
    ROI_   = UnitRange.( out_0, out_1 ); 
    WHD    = 2 .* rad

    @inbounds for c in CartesianIndices( ROI_ )
        output[ c ] = op( output[c], f * integralSums_unsafe( intArr, Tuple(c) .- rad, Tuple(c) .+ rad ) )
    end

    return nothing
end

"""
    Similar as the above function, but uses 'safe' integral sums
    and includes a paramter that controls how to combine local 
    sums and 'N'. 
"""

# r = (10,10); inp = rand( Float32, 100 + 2*r[1], 100 + 2*r[2] ); 
# @btime out = localsums_unsafe( $inp, $r )  
#   47.300 μs (4 allocations: 113.59 KiB)
function localsums_N_f_op( input::AbstractArray{T,N}, 
                           rad::Dims{N},
                           Nop::Function=*,
                           f::T=T(1),
                           op::Function=+
                              ) where {T<:AbstractFloat,N}
    output = zeros( T, size(input));
    intA   = integralArray( input ); 
    localsums_N_f_op!( intA, output, rad, Nop, f, op );
    return output;
end

# 2D
function localsums_N_f_op!( intArr::AbstractArray{T,2}, 
                            output::AbstractArray{T,2}, 
                            rad::Dims{2},
                            Nop::Function=*, 
                            f::T=T(1),
                            op::Function=+ ) where {T<:AbstractFloat}

    @inbounds for c in CartesianIndices( output )
        output[ c ] = op( output[c], f * integralSumN( intArr, Tuple(c) .- rad, Tuple(c) .+ rad, Nop ) )
    end
    return nothing
end

# 3D
function localsums_N_f_op!( intArr::AbstractArray{T,3}, 
                            output::AbstractArray{T,3}, 
                            rad::Dims{3},
                            Nop::Function=*,
                            f::T=T(1),
                            op::Function=+ )  where {T<:AbstractFloat}

    @inbounds for c in CartesianIndices( output )
        output[ c ] = op( output[c], f * integralSumsN( intArr, Tuple(c) .- rad, Tuple(c) .+ rad, Nop ) )
    end

    return nothing
end

"""
   Sometimes it is desired to divide by N or multiply by N. THe
   functions below enable that.
""" 

function localN( img, rad )
    output = zeros( eltype(img), size(img) );
    for c in CartesianIndices( img )
        TL = clipmin.( Tuple(c) .- rad, 1 )
        BR = clipmax.( Tuple(c) .+ rad, size(img) ); 
        output[c] = prod( BR .- TL .+ 1 ); 
    end
    return output
end

function localN_op!( output::AbstractArray{T,N}, 
                     rad::Dims{N}, 
                     op1::Function, 
                     op2::Function ) where {T<:AbstractFloat,N}
    for c in CartesianIndices( output )
        TL = clipmin.( Tuple(c) .- rad, 1 )
        BR = clipmax.( Tuple(c) .+ rad, size(output) ); 
        output[c] = op2( output[c],  op1( T( prod( BR .- TL .+ 1 ) ) )  ); 
    end
    return nothing
end
