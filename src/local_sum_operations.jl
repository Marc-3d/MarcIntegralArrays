##### BOUND-SAFE LOCAL SUMS

function localsums( input::AbstractArray{T,N}, 
                    rad::Dims{N} 
                   ) where {T<:AbstractFloat,N}

    output = zeros( T, size(input));
    intA   = integralArray( input ); 
    localsums!( intA, output, rad );
    return output;
end

function localsums_N_f_op( input::AbstractArray{T,N}, 
                           rad::Dims{N},
                           Nop::Function=(sum::T,n::T)->(sum/n),
                           f::T=T(1),
                           op::Function=(out::T,res::T)->(out+res)
                          ) where {T<:AbstractFloat,N}

    output = zeros( T, size(input));
    intA   = integralArray( input ); 
    localsums_N_f_op!( intA, output, rad, Nop, f, op );
    return output;
end

#= 
NOTE: These Benchmark results should be proportional to the number of elements in the input

using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; T = Float32; IA = rand( T, N+1 ); output = zeros( T, N ); rad = (5,);
@btime MarcIntegralArrays.localsums!( $IA, $output, $rad )
  ~ 250 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; T = Float32; IA = rand( T, Ny+1, Nx+1 ); output = zeros( T, Ny, Nx ); rad = (5,5);
@btime MarcIntegralArrays.localsums!( $IA, $output, $rad )
  ~ 70.600 μs (AMD EPYC 7453)

# 3D
Ny, Nx, Nz = 20, 20, 20; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); output = zeros( T, Ny, Nx, Nz ); rad = (5,5,5);
@btime MarcIntegralArrays.localsums!( $IA, $output, $rad )
   ~ 88 μs (AMD EPYC 7453)
=#
""" 
    Computes in-place Local sums each each coordinate, p, of the input data. 
    These local sums are computed from the local rectangular ROI around each
    coordiante, p, given by `UnitRange.( p.-rad, p.+rad )`.
"""
function localsums!( intArr::AbstractArray{T,N},
                     output::AbstractArray{T,N},
	                 rad::Dims{N} ) where {T,N}

    @inbounds for c in CartesianIndices( output )
        output[ c ] = integralSum( intArr, Tuple(c) .- rad, Tuple(c) .+ rad )
    end 
    return nothing
end

# in-place Local sums in a "ring ROI" (by sustracting a small rectangle from a big rectangle).
function localsums!( intArr::AbstractArray{T,N},
                     output::AbstractArray{T,N},
	                 rad_in::Dims{N},
                     rad_out::Dims{N} ) where {T,N}

    @inbounds for c in CartesianIndices( output )
        output[ c ] = integralSum( intArr, Tuple(c) .- rad_out, Tuple(c) .+ rad_out ) - integralSum( intArr, Tuple(c) .- rad_in, Tuple(c) .+ rad_in )
    end 
    return nothing
end

# IA = rand( Float32, 50, 50 ); out = zeros( Float32, 49, 49 ); rad = (3,3); nop=*; f=Float32(1); op=+;
# @btime MiA.localsums_N_f_op!( $IA, $out, $rad, $nop, $f, $op )
#   122.300 μs (7205 allocations: 187.62 KiB) (BOBA server) (vs 123.000)
#
# IA = rand( Float32, 20, 20, 20 ); out = zeros( Float32, 19, 19, 19 ); rad = (3,3,3); nop=*; f=Float32(1); op=+;
# @btime MiA.localsums_N_f_op!( $IA, $out, $rad, $nop, $f, $op )
#   1.673 ms (20579 allocations: 535.91 KiB) (BOBA server) (vs 1.674.00)
#
# IA = rand( Float32, 20, 20, 20 ); out = zeros( Float32, 19, 19, 19 ); rad = (3,3,3); nop=(x::Float32,y::Float32)->(x*y); f=Float32(1); op=(x::Float32,y::Float32)->(x+y);
# @btime MiA.localsums_N_f_op!( $IA, $out, $rad, $nop, $f, $op )
#   369.800 μs (20579 allocations: 535.91 KiB) (BOBA server) 
function localsums_N_f_op!( intArr::AbstractArray{T,N}, 
                            output::AbstractArray{T,N}, 
                            rad::Dims{N},
                            Nop::Function=(sum::T,n::T)->(sum/n), 
                            f::T=T(1),
                            op::Function=(out::T,res::T)->(out+res) 
                          ) where {T<:AbstractFloat,N}

    @inbounds for c in CartesianIndices( output )
        output[ c ] = op( output[c], f * integralSumN( intArr, Tuple(c) .- rad, Tuple(c) .+ rad, Nop ) )
    end
    return nothing
end

function localsums_N_f_op!( intArr::AbstractArray{T,N}, 
                            output::AbstractArray{T,N}, 
                            rad::Dims{N};
                            f::T=T(1),
                            Nop::Function=(sum::T,n::T)->(sum/n), 
                            op::Function=(out::T,res::T)->(out+res) 
                          ) where {T<:AbstractFloat,N}

    localsums_N_f_op!( intArr, output, rad, Nop, f, op )
    return nothing
end

function localsums_N_f_op!( intArr::AbstractArray{T,N}, 
                            output::AbstractArray{T,N}, 
                            rad1::Dims{N},
                            rad2::Dims{N},
                            Nop::Function=(sum::T,n::T)->(sum/n), 
                            f::T=T(1),
                            op::Function=(out::T,res::T)->(out+res) 
                          ) where {T<:AbstractFloat,N}

    @inbounds for c in CartesianIndices( output )
        output[ c ] = op( output[c], f * integralSumN( intArr, Tuple(c) .- rad1, Tuple(c) .+ rad1, Tuple(c) .- rad2, Tuple(c) .+ rad2, Nop ) )
    end
    return nothing
end

function localsums_N_f_op!( intArr::AbstractArray{T,N}, 
                            output::AbstractArray{T,N}, 
                            rad1::Dims{N},
                            rad2::Dims{N};
                            f::T=T(1),
                            Nop::Function=(sum::T,n::T)->(sum/n), 
                            op::Function=(out::T,res::T)->(out+res) 
                          ) where {T<:AbstractFloat,N}

    localsums_N_f_op!( intArr, output, rad1, rad2, Nop, f, op )
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

"""
    Similar as the above function, but uses 'safe' integral sums
    and includes a paramter that controls how to combine local 
    sums and 'N'. 
"""



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
                     Nop::Function, 
                     Oop::Function 
                   ) where {T<:AbstractFloat,N}

    for c in CartesianIndices( output )
        TLF = clipmin.( Tuple(c) .- rad, 1 )
        BRB = clipmax.( Tuple(c) .+ rad, size(output) ); 
        output[c] = Oop( output[c],  Nop( T( prod( BRB .- TLF .+ 1 ) ) )  ); 
    end
    return nothing
end

function localN_op!( output::AbstractArray{T,N}, 
                     rad::Dims{N};
                     Nop::Function=(n::T)->(n), 
                     op::Function=(out::T,res::T)->(out*res)
                   ) where {T<:AbstractFloat,N}
    localN_op!( output, rad, Nop, op )
    return nothing
end

function localN_op!( output::AbstractArray{T,N}, 
                     rad1::Dims{N},
                     rad2::Dims{N}, 
                     Nop::Function=(n::T)->(n), 
                     op::Function=(out::T,res::T)->(out*res)
                   ) where {T<:AbstractFloat,N}

    for c in CartesianIndices( output )
        TLF1  = clipmin.( Tuple(c) .- rad1, 1 )
        BRB1  = clipmax.( Tuple(c) .+ rad1, size(output) ); 
        TLF2  = clipmin.( Tuple(c) .- rad2, 1 )
        BRB2  = clipmax.( Tuple(c) .+ rad2, size(output) ); 
        Nring = T( prod( BRB2 .- TLF2 .+ 1 ) - prod( BRB1 .- TLF1 .+ 1 ) )
        output[c] = op( output[c],  Nop( Nring )  ); 
    end
    return nothing
end

function localN_op!( output::AbstractArray{T,N}, 
                     rad1::Dims{N},
                     rad2::Dims{N};
                     Nop::Function=(n::T)->(n), 
                     op::Function=(out::T,res::T)->(out*res)
                   ) where {T<:AbstractFloat,N}
    localN_op!( output, rad1, rad2, Nop, op )
    return nothing
end
