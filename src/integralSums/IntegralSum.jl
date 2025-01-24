#=
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; IA = rand( Float32, N+1 ); L = rand( 1:N ); R = rand(L+1:N); f = T(3);
@btime MarcIntegralArrays.integralSum( $IA, $L, $R, $f )
  ~ 8 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; IA = rand( Float32, Ny+1, Nx+1 ); TL = rand.( (1:Ny,1:Nx) ); BR = rand.((TL[1]+1:Ny,TL[2]+1:Nx)); f = T(3);
@btime MarcIntegralArrays.integralSum( $IA, $TL, $BR, $f )
  ~ 9-10 ns (AMD EPYC 7453)

# 3D
Ny, Nx, Nz = 30, 30, 30; IA = rand( Float32, Ny+1, Nx+1, Nz+1 ); TLF = rand.( (1:Ny,1:Nx,1:Nz) ); BRB = rand.((TLF[1]+1:Ny,TLF[2]+1:Nx,TLF[3]+1:Nz)); f = T(3);
@btime MarcIntegralArrays.integralSum( $IA, $TLF, $BRB, $f )
  ~ 14 ns (AMD EPYC 7453)
=#
"""
  Computes the sum of intensities within a rectangular ROI defined by its top-left[-front] (TL[F]) and bottom-right[-back] (BR[B]) corners. This is a 'safe' function (in contrast to `integralSum_unsafe`), which crops the input ROI to make sure it doesn't go out-of-bounds. By using the inlined functions above, the 'safe' code is suprisingly on-par with the unsafe code in terms of performance. Therefore, you probably should stick to the 'safe' function, unless you are very desperate for performance.
"""
function integralSum( 
  intArr::AbstractArray{T,N}, 
  TL::Dims{N},
  BR::Dims{N},
  f::T=T(1)
) where {
  T,
  N
}
  return integralSum_unsafe( 
    intArr, 
    minmax( TL, 1, size(intArr).-1 ), 
    minmax( BR, 1, size(intArr).-1 ),
    f
  )
end

# 1D: using Int instead of of Dims{1}
function integralSum( 
  intArr::AbstractArray{T,1}, 
  L::Int,
  R::Int,
  f::T=T(1)
) where {
  T
}
  return integralSum_unsafe( 
    intArr, 
    minmax( L, 1, size(intArr,1)-1 ), 
    minmax( R, 1, size(intArr,1)-1 ),
    f
  )
end

############

#=
using BenchmarkTools, MarcIntegralArrays

N = 100; IA = rand( Float32, N+1 ); L = rand( 1:N ); R = rand(L+1:N);
@btime MarcIntegralArrays.integralSum_unsafe( $IA, $L, $R )
  ~ 8 ns (AMD EPYC 7453)

Ny, Nx = 100, 100; IA = rand( Float32, Ny+1, Nx+1 ); TL = rand.( (1:Ny,1:Nx) ); BR = rand.((TL[1]+1:Ny,TL[2]+1:Nx));
@btime MarcIntegralArrays.integralSum_unsafe( $IA, $TL, $BR )
  ~ 8 ns (AMD EPYC 7453)

Ny, Nx, Nz = 30, 30, 30; IA = rand( Float32, Ny+1, Nx+1, Nz+1 ); TLF = rand.( (1:Ny,1:Nx,1:Nz) ); BRB = rand.((TLF[1]+1:Ny,TLF[2]+1:Nx,TLF[3]+1:Nz)); f = T(3)
@btime MarcIntegralArrays.integralSum_unsafe( $IA, $TLF, $BRB, $f )
  ~ 10 ns (AMD EPYC 7453)
=#

""" 1D integrals sums unsafe """
function integralSum_unsafe( 
  intArr::AbstractArray{T,1}, 
  L::Dims{1}, 
  R::Dims{1},
  f::T=T(1)
) where {
  T<:AbstractFloat
}
  return integralSum_unsafe( 
    intArr, 
    L[1], 
    R[1],
    f 
  )
end

function integralSum_unsafe( 
  intArr::AbstractArray{T,1}, 
  L::Int, 
  R::Int,
  f::T=T(1)
)::T where {
  T<:AbstractFloat
}
  @inbounds begin  
      return f * ( intArr[ R+1 ]  - intArr[ L ] )
  end
end

""" 2D integral sums unsafe """
function integralSum_unsafe( 
  intArr::AbstractArray{T,2}, 
  TL::Dims{2}, 
  BR::Dims{2},
  f::T=T(1)
)::T where {
  T<:AbstractFloat
}
  @inbounds begin  
      return f * ( intArr[TL[1],TL[2]] - intArr[BR[1]+1,TL[2]] - intArr[TL[1],BR[2]+1] + intArr[BR[1]+1,BR[2]+1] )
  end
end

""" 3D integral sums unsafe """
function integralSum_unsafe( 
  intArr::AbstractArray{T,3}, 
  TLF::Dims{3}, 
  BRB::Dims{3},
  f::T=T(1)
)::T where {
  T<:AbstractFloat
}
  @inbounds begin 
    return f * ( intArr[BRB[1]+1,BRB[2]+1,BRB[3]+1] - intArr[TLF[1],TLF[2],TLF[3]] - intArr[TLF[1],BRB[2]+1,BRB[3]+1] - intArr[BRB[1]+1,TLF[2],BRB[3]+1] - intArr[BRB[1]+1,BRB[2]+1,TLF[3]] + intArr[BRB[1]+1,TLF[2],TLF[3]] + intArr[TLF[1],BRB[2]+1,TLF[3]] + intArr[TLF[1],TLF[2],BRB[3]+1] )
  end    
end