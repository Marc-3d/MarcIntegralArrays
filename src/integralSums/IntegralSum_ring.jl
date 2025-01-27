#=
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; IA = rand( Float32, N+1 ); L1 = rand( 1:N ); R1 = rand(L1+1:N); L2 = rand( 1:N ); R2 = rand(L2+1:N);
@btime MarcIntegralArrays.integralSum( $IA, $L1, $R1, $L2, $R2 )
  ~ 13 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; IA = rand( Float32, Ny+1, Nx+1 ); TL1 = rand.( (1:Ny,1:Nx) ); BR1 = rand.((TL1[1]+1:Ny,TL1[2]+1:Nx)); TL2 = rand.( (1:Ny,1:Nx) ); BR2 = rand.((TL2[1]+1:Ny,TL2[2]+1:Nx));
@btime MarcIntegralArrays.integralSum( $IA, $TL1, $BR1, $TL2, $BR2 )
  ~ 20 ns (AMD EPYC 7453)

# 3D
Ny, Nx, Nz = 30, 30, 30; IA = rand( Float32, Ny+1, Nx+1, Nz+1 ); TLF1 = rand.( (1:Ny,1:Nx,1:Nz) ); BRB1 = rand.((TLF1[1]+1:Ny,TLF1[2]+1:Nx,TLF1[3]+1:Nz)); TLF2 = rand.( (1:Ny,1:Nx,1:Nz) ); BRB2 = rand.((TLF2[1]+1:Ny,TLF2[2]+1:Nx,TLF2[3]+1:Nz));
@btime MarcIntegralArrays.integralSum( $IA, $TLF1, $BRB1, $TLF2, $BRB2 )
  ~ 30 ns (AMD EPYC 7453)
=#

"""
    This function accepts two radii, and computes the sum in a ring by substracting the sums in small ROI to the sums in a bigger ROI. The first set of ROI coordinates (TL1 and BR1) define the smaller ROI, and the TL2 and BR2 define the larger ROI. It basically involves computing two integral sums, and it is about twice as slow as computing a single integral sum.
"""
function integralSum( 
  intArr::AbstractArray{T,N}, 
  TL1::Dims{N},
  BR1::Dims{N},
  TL2::Dims{N},
  BR2::Dims{N},
  f::T=T(1)
) where {
  T,
  N
}
    return integralSum_unsafe( 
        intArr, 
        minmax( TL2, 1, size(intArr).-1 ), 
        minmax( BR2, 1, size(intArr).-1 ),
        f
    ) - integralSum_unsafe(
        intArr, 
        minmax( TL1, 1, size(intArr).-1 ), 
        minmax( BR1, 1, size(intArr).-1 ),
        f  
    )
end

# 1D: using Int instead of of Dims{1}
function integralSum(
  intArr::AbstractArray{T,1}, 
  L1::Int,
  R1::Int,
  L2::Int,
  R2::Int,
  f::T=T(1) 
) where {
  T
}
    return integralSum_unsafe( 
        intArr, 
        minmax( L2, 1, size(intArr,1)-1 ), 
        minmax( R2, 1, size(intArr,1)-1 ),
        f
    ) - integralSum_unsafe(
        intArr,
        minmax( L1, 1, size(intArr,1)-1 ), 
        minmax( R1, 1, size(intArr,1)-1 ),
        f
    )
end
