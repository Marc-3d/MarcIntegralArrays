# NOTE: All Benchmarks below are proportional to the number of elements in the input. In other words, it is expected that computing the local sum around each voxel in a 100x100x100 volume (1 million voxels) will take much longer (100 times longer) than computing the local sum around each pixel in a 100x100 image (10 thousand pixels)
# TODO: I removed some functions from this file (localN_op! and localsums_N_f_op). Make sure these functions aren't used anywere else

##### BOUND-SAFE LOCAL SUMS

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
    NN     = prod( 2 .* rad .+ 1 ) 
    out_v  = view( output, :, :, 1 )
    localsums!( out_v, intA.arr, rad );

    for i in 2:3
        c = channels[i]
        integralArray!( intA, input, channel=c ); 
        out_v = view( output, :, :, c )
        localsums!( out_v, intA.arr, rad )
    end
    output ./= NN

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
    NN     = prod( 2 .* rad .+ 1 ) 
    out_v  = view( output, :, :, 1 )
    localsums!( out_v, intA.arr, rad );

    for i in 2:3
        c = channels[i]
        integralArray!( intA, input, channel=c ); 
        out_v  = view( output, :, :, c )
        localsums!( out_v, intA.arr, rad )
    end
    # TOOD: hard-corded for 2D at the moment
    for c in CartesianIndices( input )
        output[Tuple(c)...,:] ./= prod( 2 .* (rad[c],rad[c]) .+ 1 )
    end

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
    Computes in-place local sums for each coordinate, p, of the input data. These local sums are computed from the local rectangular ROI around each coordiante, p, given by `UnitRange.( p.-rad, p.+rad )`.
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
            # TODO: hard-coded for 2D at the moment. the "rad" tuple should have N elements.
            Tuple(c) .- (rad[c],rad[c]), 
            Tuple(c) .+ (rad[c],rad[c]),
            f 
        )
        output[ c ] = op( output[c], tmp ); 
    end 
    return nothing
end



#= 
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; T = Float32; IA = rand( T, N+1 ); output = zeros( T, N ); rad_in = (5,); rad_out=(10,); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localsums!( $IA, $output, $rad_in, $rad_out, f=$f, op=$op )
  ~ 500 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; T = Float32; IA = rand( T, Ny+1, Nx+1 ); output = zeros( T, Ny, Nx ); rad_in = (5,5); rad_out=(10,10); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localsums!( $IA, $output, $rad_in, $rad_out, f=$f, op=$op )
  ~ 140.600 μs (AMD EPYC 7453)

# 3D
Ny, Nx, Nz = 20, 20, 20; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); output = zeros( T, Ny, Nx, Nz ); rad_in = (5,5,5); rad_out = (10,10,10); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localsums!( $IA, $output, $rad_in, $rad_out, f=$f, op=$op )
   ~ 170 μs (AMD EPYC 7453)
=#
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

#= 
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; T = Float32; IA = rand( T, N+1 ); output = zeros( T, N ); rad = (5,); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localsumNs!( $IA, $output, $rad, f=$f, op=$op )
  ~ 250 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; T = Float32; IA = rand( T, Ny+1, Nx+1 ); output = zeros( T, Ny, Nx ); rad = (5,5); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localsumNs!( $IA, $output, $rad, f=$f, op=$op )
  ~ 260.600 μs (AMD EPYC 7453)

# 3D
Ny, Nx, Nz = 20, 20, 20; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); output = zeros( T, Ny, Nx, Nz ); rad = (5,5,5); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localsumNs!( $IA, $output, $rad, f=$f, op=$op )
   ~ 300 μs (AMD EPYC 7453)
=#
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

#= 
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; T = Float32; IA = rand( T, N+1 ); output = zeros( T, N ); rad = (5,); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localAvgs!( $IA, $output, $rad, f=$f, op=$op )
  ~ 250 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; T = Float32; IA = rand( T, Ny+1, Nx+1 ); output = zeros( T, Ny, Nx ); rad = (5,5); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localAvgs!( $IA, $output, $rad, f=$f, op=$op )
  ~ 260.600 μs (AMD EPYC 7453)

# 3D
Ny, Nx, Nz = 20, 20, 20; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); output = zeros( T, Ny, Nx, Nz ); rad = (5,5,5); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
@btime MarcIntegralArrays.localAvgs!( $IA, $output, $rad, f=$f, op=$op )
   ~ 300 μs (AMD EPYC 7453)
=#
""" 
    Computes in-place local average ( local sums / number of pixels ) in the rectangular local ROI.
"""
function localAvgs!( 
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

function localAvgs( 
    input::AbstractArray{T,N},
    rad::Dims{N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T,
    N
}
    intA = IntegralArray( input )
    output = zeros( T, size(input) )
    localAvgs!( output, intA.arr, rad, f=f, op=op )
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

function localAvgs( 
    input::AbstractArray{T,N},
    rad1::Dims{N},
    rad2::Dims{N};
    f::T=T(1),
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T,
    N
}
    intA = integralArray( input )
    output = zeros( T, size(input) )
    localAvgs!( output, intA, rad1, rad2, f=f, op=op )
    return output
end

"""
   Sometimes it is desired to divide by N or multiply by N. The functions below enable that.
""" 
function localN( 
    img::AbstractArray{T,N}, 
    rad::Dims{N}
) where {
    T,
    N
}
    output = zeros( eltype(img), size(img) );
    localN!( output, rad )
    return output
end

#= 
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; T = Float32; output = zeros( T, N ); rad = (5,); f = Float32(1); op = (o::Float32,n::Float32)->(o/n)::Float32;
@btime MarcIntegralArrays.localN!( $output, $rad, f=$f, op=$op )
  ~ 100 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; T = Float32; output = zeros( T, Ny, Nx ); rad = (5,5); f = Float32(1); op = (o::Float32,n::Float32)->(o/n)::Float32;
@btime MarcIntegralArrays.localN!( $output, $rad, f=$f, op=$op )
  ~ 12.0 μs (AMD EPYC 7453)

# 3D
Ny, Nx, Nz = 20, 20, 20; T = Float32; output = zeros( T, Ny, Nx, Nz ); rad = (5,5,5); f = Float32(1); op = (o::Float32,n::Float32)->(o/n)::Float32;
@btime MarcIntegralArrays.localN!( $output, $rad, f=$f, op=$op )
   ~ 53 μs (AMD EPYC 7453)
=#
function localN!(
    output::AbstractArray{T,N}, 
    rad::Dims{N};
    f::T=T(1),
    op=(o::T,n)->(o/n)
) where {
    T,
    N
}
    for c in CartesianIndices( output )
        TL = clipmin.( Tuple(c) .- rad, 1 )
        BR = clipmax.( Tuple(c) .+ rad, size(output) ); 
        N_ = f * T( prod( BR .- TL .+ 1 ) )
        output[c] = op( output[c], N_ ); 
    end
    return nothing
end

function localN!(
    output::AbstractArray{T,N}, 
    rad_in::Dims{N},
    rad_out::Dims{N};
    f::T=T(1),
    op=(o::T,n)->(o/n)
) where {
    T,
    N
}
    for c in CartesianIndices( output )
        TL1  = clipmin.( Tuple(c) .- rad_in, 1 )
        BR1  = clipmax.( Tuple(c) .+ rad_in, size(output) ); 
        N_in = f * T( prod( BR1 .- TL1 .+ 1 ) )

        TL2  = clipmin.( Tuple(c) .- rad_out, 1 )
        BR2  = clipmax.( Tuple(c) .+ rad_out, size(output) ); 
        Nout = f * T( prod( BR2 .- TL2 .+ 1 ) )

        output[c] = op( output[c], Nout - N_in ); 
    end
    return nothing
end

########## UNSAFE FUNCTIONS

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

# IA = rand( Float32, 50, 50 ); out = zeros( Float32, 49, 49 ); rad = (3,3);
# @btime localsums_unsafe_ND!( $IA, $output, $rad )
#  1.700 μs (1 allocation: 32 bytes) (BOBA server) (vs 1.750)
#
# IA = rand( Float32, 50, 50, 50 ); out = zeros( Float32, 49, 49, 49 ); rad = (3,3, 3);
# @btime localsums_unsafe!( $IA, $output, $rad )
#  303.700 μs (1 allocation: 32 bytes) (BOBA server) (vs 303.600)
function localsums_unsafe!( intArr::AbstractArray{T,N}, 
                            output::AbstractArray{T,N}, 
                            rad::Dims{N} ) where {T,N}

    out_sz = size( output )
    out_0  = rad .+ 1
    out_1  = out_sz .- rad
    ROI_   = UnitRange.( out_0, out_1 ); 

    @inbounds for c in CartesianIndices( ROI_ )
        output[ c ] = integralSum_unsafe( intArr, Tuple(c) .- rad, Tuple(c) .+ rad )
    end
    return nothing
end

#=
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
=#

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
    out_0  = rad .+ 1
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
    out_0  = rad .+ 1
    out_1  = out_sz  .- rad
    ROI_   = UnitRange.( out_0, out_1 ); 
    WHD    = 2 .* rad

    @inbounds for c in CartesianIndices( ROI_ )
        output[ c ] = op( output[c], f * integralSums_unsafe( intArr, Tuple(c) .- rad, Tuple(c) .+ rad ) )
    end

    return nothing
end
