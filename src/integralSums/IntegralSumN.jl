#=
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; T = Float32; IA = rand( T, N+1 ); L = rand( 1:N ); R = rand(L+1:N); f = T(2);
@btime MarcIntegralArrays.integralSumN( $IA, $L, $R, $f )
@btime MarcIntegralArrays.integralAvg( $IA, $L, $R, $f )
  ~ 8 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; T = Float32; IA = rand( T, Ny+1, Nx+1 ); TL = rand.( (1:Ny,1:Nx) ); BR = rand.((TL[1]+1:Ny,TL[2]+1:Nx)); f = T(2);
@btime MarcIntegralArrays.integralSumN( $IA, $TL, $BR, $f )
@btime MarcIntegralArrays.integralAvg( $IA, $TL, $BR, $f )
  ~ 14 ns (AMD EPYC 7453) 

# 3D
Ny, Nx, Nz = 30, 30, 30; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); TLF = rand.( (1:Ny,1:Nx,1:Nz) ); BRB = rand.((TLF[1]+1:Ny,TLF[2]+1:Nx,TLF[3]+1:Nz)); f = T(2); 
@btime MarcIntegralArrays.integralSumN( $IA, $TLF, $BRB, $f )
@btime MarcIntegralArrays.integralAvg( $IA, $TLF, $BRB, $f )
  ~ 22 ns (AMD EPYC 7453)
=#
"""
    Computes integral sum multiplied by the number of elements in the ROI.
"""
function integralSumN( 
    intArr::AbstractArray{T,N},
    TL::Dims{N},
    BR::Dims{N},
    f::T=T(1)
)::T where {
    T,
    N
}
  return integralSumN_unsafe( 
    intArr, 
    minmax( TL, 1, size(intArr).-1 ), 
    minmax( BR, 1, size(intArr).-1 ),
    f
  );
end

function integralSumN( 
    intArr::AbstractArray{T,1}, 
    L::Int,
    R::Int, 
    f::T=T(1)
)::T where {
    T
}
  return integralSumN_unsafe( 
    intArr, 
    minmax( L, 1, size(intArr,1)-1 ), 
    minmax( R, 1, size(intArr,1)-1 ),
    f
  )
end

"""
    Computes integral sum divided by the number of elements in the ROI
"""
function integralAvg(
    intArr::AbstractArray{T,N}, 
    TL::Dims{N}, 
    BR::Dims{N},
    f::T=T(1)
) where {
    T,
    N
}
    return integralAvg_unsafe( 
        intArr, 
        minmax( TL, 1, size(intArr).-1 ), 
        minmax( BR, 1, size(intArr).-1 ),
        f
    )
end

function integralAvg(
    intArr::AbstractArray{T,1}, 
    L::Int, 
    R::Int,
    f::T=T(1)
)::T where {
    T
}
    return integralAvg_unsafe( 
        intArr, 
        minmax( L, 1, size(intArr,1)-1 ), 
        minmax( R, 1, size(intArr,1)-1 ),
        f
    )
end

####################################################### 

#=
using BenchmarkTools, MarcIntegralArrays

N = 100; T=Float32; IA = rand( T, N+1 ); L = rand( 1:N ); R = rand(L+1:N);
@btime MarcIntegralArrays.integralSumN_unsafe( $IA, $L, $R )
@btime MarcIntegralArrays.integralAvg_unsafe( $IA, $L, $R )
  ~ 8 ns (AMD EPYC 7453)

Ny, Nx = 100, 100; IA = rand( Float32, Ny+1, Nx+1 ); TL = rand.( (1:Ny,1:Nx) ); BR = rand.((TL[1]+1:Ny,TL[2]+1:Nx)); 
@btime MarcIntegralArrays.integralSumN_unsafe( $IA, $TL, $BR )
@btime MarcIntegralArrays.integralAvg_unsafe( $IA, $TL, $BR )
  ~ 8 ns (AMD EPYC 7453)

Ny, Nx, Nz = 30, 30, 30; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); TLF = rand.( (1:Ny,1:Nx,1:Nz) ); BRB = rand.((TLF[1]+1:Ny,TLF[2]+1:Nx,TLF[3]+1:Nz)); 
@btime MarcIntegralArrays.integralSumN_unsafe( $IA, $TLF, $BRB )
@btime MarcIntegralArrays.integralAvg_unsafe( $IA, $TLF, $BRB )
  ~ 13 ns (AMD EPYC 7453)
=#

""" 1D """
function integralSumN_unsafe( 
    intArr::AbstractArray{T,1}, 
    L::Int, 
    R::Int, 
    f::T=T(1)
)::T where {
    T<:AbstractFloat
}
    @inbounds begin  
        return f * ( intArr[ R+1 ] - intArr[ L ] ) * T( R - L + 1 ) 
    end
end

""" 1D, 2D, 3D """
function integralSumN_unsafe( 
    intArr::AbstractArray{T,N}, 
    TL::Dims{N}, 
    BR::Dims{N},
    f::T=T(1)
)::T where {
    T<:AbstractFloat, 
    N
}
    @inbounds begin  
        NN = prod( BR .- TL .+ 1 ) 
        return integralSum_unsafe( intArr, TL, BR, f ) * T( NN )
    end
end

####################################################### 

""" 1D """
function integralAvg_unsafe( 
    intArr::AbstractArray{T,1}, 
    L::Int, 
    R::Int,
    f::T=T(1) 
)::T where {
    T<:AbstractFloat
}
    @inbounds begin  
        return f * ( intArr[ R+1 ] - intArr[ L ] ) / T( R - L + 1 ) 
    end
end

""" 1D, 2D, 3D """
function integralAvg_unsafe( 
    intArr::AbstractArray{T,N}, 
    TL::Dims{N}, 
    BR::Dims{N},
    f::T=T(1)
)::T where {
    T<:AbstractFloat, 
    N
}
    @inbounds begin  
        NN = prod( BR .- TL .+ 1 ) 
        return integralSum_unsafe( intArr, TL, BR, f ) / T( NN )
    end
end