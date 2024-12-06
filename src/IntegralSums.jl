@inline  minmax( a, min=1, max=10 ) = a + (min-a)*(a<min) + (a>max)*(max-a)
@inline clipmin( a, min=1 ) = a + (min-a)*(a<min)
@inline clipmax( a, max=1 ) = a + (a>max)*(max-a)

#=
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; IA = rand( Float32, N+1 ); L = rand( 1:N ); R = rand(L+1:N);
@btime MarcIntegralArrays.integralSum( $IA, $L, $R )
  ~ 8 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; IA = rand( Float32, Ny+1, Nx+1 ); TL = rand.( (1:Ny,1:Nx) ); BR = rand.((TL[1]+1:Ny,TL[2]+1:Nx));
@btime MarcIntegralArrays.integralSum( $IA, $TL, $BR )
  ~ 9-10 ns (AMD EPYC 7453)

# 3D
Ny, Nx, Nz = 30, 30, 30; IA = rand( Float32, Ny+1, Nx+1, Nz+1 ); TLF = rand.( (1:Ny,1:Nx,1:Nz) ); BRB = rand.((TLF[1]+1:Ny,TLF[2]+1:Nx,TLF[3]+1:Nz));
@btime MarcIntegralArrays.integralSum( $IA, $TLF, $BRB )
  ~ 14 ns (AMD EPYC 7453)
=#

"""
   Computes the sum of intensities within a rectangular ROI defined by its 
  top-left[-front] (TL[F]) and bottom-right[-back] (BR[B]) corners. This is a
  'safe' function (in contrast to `integralSum_unsafe`), which crops the input
  ROI to make sure it doesn't go out-of-bounds. By using the inlined functions 
  above, the 'safe' code is suprisingly on-par with the unsafe code in terms of
  performance. Therefore, you probably should stick to the 'safe' function, unless
  you are very desperate for performance.
"""
function integralSum( intArr::AbstractArray{T,N}, 
                      TL::Dims{N},
                      BR::Dims{N} 
                    ) where {T,N}

    return integralSum_unsafe( intArr, 
                               minmax.( TL, 1, size(intArr).-1 ), 
                               minmax.( BR, 1, size(intArr).-1 ) );
end

function integralSum( intArr::AbstractArray{T,1}, 
                      L::Int,
                      R::Int 
                    ) where {T}

    return integralSum_unsafe( intArr, 
                               minmax( L, 1, size(intArr,1)-1 ), 
                               minmax( R, 1, size(intArr,1)-1 ) );
end

#=
using BenchmarkTools, MarcIntegralArrays

#1D
N = 100; T = Float32; IA = rand( T, N+1 ); L = rand( 1:N ); R = rand(L+1:N); op=(sum::T,N::T)->(sum/N)
@btime MarcIntegralArrays.integralSumN( $IA, $L, $R, $op )
  ~ 7 ns (AMD EPYC 7453)

# 2D
Ny, Nx = 100, 100; T = Float32; IA = rand( T, Ny+1, Nx+1 ); TL = rand.( (1:Ny,1:Nx) ); BR = rand.((TL[1]+1:Ny,TL[2]+1:Nx)); op=(sum::T,N::T)->(sum/N)
@btime MarcIntegralArrays.integralSumN( $IA, $TL, $BR, $op )
  ~ 50 ns (AMD EPYC 7453) TODO: FIX ALLOCATIONS IN DENOMINATOR

# 3D
Ny, Nx, Nz = 30, 30, 30; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); TLF = rand.( (1:Ny,1:Nx,1:Nz) ); BRB = rand.((TLF[1]+1:Ny,TLF[2]+1:Nx,TLF[3]+1:Nz)); op=(sum::T,N::T)->(sum/N)
@btime MarcIntegralArrays.integralSumN( $IA, $TLF, $BRB, $op )
  ~ 80 ns (AMD EPYC 7453) TODO: FIX ALLOCATIONS IN DENOMINATOR
=#
"""
    Computes integral sum, and operates on the number of elements in the ROI.
"""
function integralSumN( intArr::AbstractArray{T,N},
                       TL::Dims{N},
                       BR::Dims{N},
                       op::Function=(sum::T,n::T)->(sum/n)
                     )::T where {T,N}

    return integralSumN_unsafe( intArr, 
                                minmax.( TL, 1, size(intArr).-1 ), 
                                minmax.( BR, 1, size(intArr).-1 ), 
                                op );
end

function integralSumN( intArr::AbstractArray{T,1}, 
                       L::Int,
                       R::Int, 
                       op::Function=(sum::T,n::T)->(sum/n) 
                     )::T where {T}

    return integralSumN_unsafe( intArr, 
                               minmax( L, 1, size(intArr,1)-1 ), 
                               minmax( R, 1, size(intArr,1)-1 ),
                               op );
end

function integralAverage( intArr::AbstractArray{T,N}, TL::Dims{N}, BR::Dims{N} ) where {T,N}
    return integralSumN( intArr, TL, BR );
end

"""
    RING
"""
function integralSum( intArr::AbstractArray{T,N}, 
                      TL1::Dims{N},
                      BR1::Dims{N},
                      TL2::Dims{N},
                      BR2::Dims{N} 
                    ) where {T,N}

    return integralSum_unsafe( intArr, 
                               minmax.( TL1, 1, size(intArr).-1 ), 
                               minmax.( BR1, 1, size(intArr).-1 ),
                               minmax.( TL2, 1, size(intArr).-1 ), 
                               minmax.( BR2, 1, size(intArr).-1 ) );
end

function integralSum( intArr::AbstractArray{T,1}, 
                      L1::Int,
                      R1::Int,
                      L2::Int,
                      R2::Int  
                    ) where {T}

    return integralSum_unsafe( intArr, 
                               minmax( L1, 1, size(intArr,1)-1 ), 
                               minmax( R1, 1, size(intArr,1)-1 ),
                               minmax( L2, 1, size(intArr,1)-1 ), 
                               minmax( R2, 1, size(intArr,1)-1 ) );
end

"""
    RING N
"""
function integralSumN( intArr::AbstractArray{T,N}, 
                      TL1::Dims{N},
                      BR1::Dims{N},
                      TL2::Dims{N},
                      BR2::Dims{N},
                      op::Function=(sum::T,n::T)->(sum/n) 
                    ) where {T,N}

    return integralSumN_unsafe( intArr, 
                               minmax.( TL1, 1, size(intArr).-1 ), 
                               minmax.( BR1, 1, size(intArr).-1 ),
                               minmax.( TL2, 1, size(intArr).-1 ), 
                               minmax.( BR2, 1, size(intArr).-1 ),
                               op );
end

function integralSumN( intArr::AbstractArray{T,1}, 
                      L1::Int,
                      R1::Int,
                      L2::Int,
                      R2::Int,
                      op::Function=(sum::T,n::T)->(sum/n) 
                    ) where {T}

    return integralSumN_unsafe( intArr, 
                               minmax( L1, 1, size(intArr,1)-1 ), 
                               minmax( R1, 1, size(intArr,1)-1 ),
                               minmax( L2, 1, size(intArr,1)-1 ), 
                               minmax( R2, 1, size(intArr,1)-1 ),
                               op );
end


####################################################### 


""" 1D """
function integralSum_unsafe( intArr::AbstractArray{T,1}, 
                             L::Dims{1}, 
                             R::Dims{1} 
                           ) where {T<:AbstractFloat}
    return integralSum_unsafe( intArr, L[1], R[1] )
end

#=
using BenchmarkTools, MarcIntegralArrays
N = 100; IA = rand( Float32, N+1 ); L = rand( 1:N ); R = rand(L+1:N);
@btime MarcIntegralArrays.integralSum_unsafe( $IA, $L, $R )
  ~ 8 ns (AMD EPYC 7453)
=#
function integralSum_unsafe( intArr::AbstractArray{T,1}, 
                             L::Int, 
                             R::Int 
                           )::T where {T<:AbstractFloat}
    @inbounds begin  
        return  intArr[ R+1 ]  - intArr[ L ] 
    end
end

""" 2D """
#=
using BenchmarkTools, MarcIntegralArrays
Ny, Nx = 100, 100; IA = rand( Float32, Ny+1, Nx+1 ); TL = rand.( (1:Ny,1:Nx) ); BR = rand.((TL[1]+1:Ny,TL[2]+1:Nx));
@btime MarcIntegralArrays.integralSum_unsafe( $IA, $TL, $BR )
  ~ 8 ns (AMD EPYC 7453)
=#
function integralSum_unsafe( intArr::AbstractArray{T,2}, 
                             TL::Dims{2}, 
                             BR::Dims{2} 
                           )::T where {T<:AbstractFloat}
    @inbounds begin  
        return intArr[TL[1], TL[2] ] - intArr[BR[1]+1, TL[2] ] - intArr[TL[1],BR[2]+1] + intArr[BR[1]+1,BR[2]+1] 
    end
end

""" 3D """
#=
using BenchmarkTools, MarcIntegralArrays
Ny, Nx, Nz = 30, 30, 30; IA = rand( Float32, Ny+1, Nx+1, Nz+1 ); TLF = rand.( (1:Ny,1:Nx,1:Nz) ); BRB = rand.((TLF[1]+1:Ny,TLF[2]+1:Nx,TLF[3]+1:Nz));
@btime MarcIntegralArrays.integralSum_unsafe( $IA, $TLF, $BRB )
  ~ 10 ns (AMD EPYC 7453)
=#
function integralSum_unsafe( intArr::AbstractArray{T,3}, 
                             TLF::Dims{3}, 
                             BRB::Dims{3} 
                           )::T where {T<:AbstractFloat}
    @inbounds begin 
        return intArr[BRB[1]+1,BRB[2]+1,BRB[3]+1] - intArr[TLF[1],TLF[2],TLF[3]] - intArr[TLF[1],BRB[2]+1,BRB[3]+1] - intArr[BRB[1]+1,TLF[2],BRB[3]+1] - intArr[BRB[1]+1,BRB[2]+1,TLF[3]] + intArr[BRB[1]+1,TLF[2],TLF[3]] + intArr[TLF[1],BRB[2]+1,TLF[3]] + intArr[TLF[1],TLF[2],BRB[3]+1]
    end    
end

# A version that is aware of "N", the size of the ROI after bound-safe cropping,
# and combined the integralSum and N according to the "two-input" function "op".

""" 1D """
function integralSumN_unsafe( intArr::AbstractArray{T,1}, 
                              L::Int, 
                              R::Int, 
                              op::Function=(sum::T,n::T)->(sum/n)
                            )::T where {T<:AbstractFloat}
    @inbounds begin  
        return op( ( intArr[ R+1 ] - intArr[ L ] ), T( R - L + 1 ) )
    end
end

#=
using BenchmarkTools, MarcIntegralArrays
N = 100; T=Float32; IA = rand( T, N+1 ); L = rand( 1:N ); R = rand(L+1:N); op=(sum::T,N::T)->(sum/N)
@btime MarcIntegralArrays.integralSumN_unsafe( $IA, $L, $R, $op )
  ~ 8 ns (AMD EPYC 7453)
=#
function integralSumN_unsafe( intArr::AbstractArray{T,1}, 
                              L::Dims{1}, 
                              R::Dims{1}, 
                              op::Function=(sum::T,n::T)->(sum/n)
                            )::T where {T<:AbstractFloat}
    return integralSumN_unsafe( intArr, L[1], R[1], op )
end

""" 2D """
#=
using BenchmarkTools, MarcIntegralArrays
Ny, Nx = 100, 100; T=Float32; IA = rand( T, Ny+1, Nx+1 ); TL = rand.( (1:Ny,1:Nx) ); BR = rand.((TL[1]+1:Ny,TL[2]+1:Nx)); op=(sum::T,N::T)->(sum/N)
@btime MarcIntegralArrays.integralSumN_unsafe( $IA, $TL, $BR, $op )
  ~ 8 ns (AMD EPYC 7453)
=#
function integralSumN_unsafe( intArr::AbstractArray{T,2}, 
                              TL::Dims{2}, 
                              BR::Dims{2}, 
                              op::Function=(sum::T,n::T)->(sum/n)
                            )::T where {T<:AbstractFloat}
    @inbounds begin  
        return op( ( intArr[TL[1],TL[2]] - intArr[BR[1]+1,TL[2]] - intArr[TL[1],BR[2]+1] + intArr[BR[1]+1,BR[2]+1] ), T( ( BR[1] - TL[1] + 1 )*( BR[2] - TL[2] + 1 ) ) )
    end
end

""" 3D """
#=
using BenchmarkTools, MarcIntegralArrays
Ny, Nx, Nz = 30, 30, 30; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); TLF = rand.( (1:Ny,1:Nx,1:Nz) ); BRB = rand.((TLF[1]+1:Ny,TLF[2]+1:Nx,TLF[3]+1:Nz)); op=(sum::T,N::T)->(sum/N)
@btime MarcIntegralArrays.integralSumN_unsafe( $IA, $TLF, $BRB, $op )
  ~ 13 ns (AMD EPYC 7453)
=#
function integralSumN_unsafe( intArr::AbstractArray{T,3}, 
                              TLF::Dims{3}, 
                              BRB::Dims{3}, 
                              op::Function=(sum::T,n::T)->(sum/n)
                            )::T where {T<:AbstractFloat}
    @inbounds begin 
        return op( ( intArr[BRB[1]+1,BRB[2]+1,BRB[3]+1] - intArr[TLF[1],TLF[2],TLF[3]] - intArr[TLF[1],BRB[2]+1,BRB[3]+1] - intArr[BRB[1]+1,TLF[2],BRB[3]+1] - intArr[BRB[1]+1,BRB[2]+1,TLF[3]] + intArr[BRB[1]+1,TLF[2],TLF[3]] + intArr[TLF[1],BRB[2]+1,TLF[3]] + intArr[TLF[1],TLF[2],BRB[3]+1] ), T( prod( BRB .- TLF .+ 1 ) ) )
    end    
end

####################################################### 

""" 1D """
function integralSum_unsafe( intArr::AbstractArray{T,1}, 
                             L1::Dims{1}, 
                             R1::Dims{1},
                             L2::Dims{1}, 
                             R2::Dims{1}
                           ) where {T<:AbstractFloat}
    return integralSum_unsafe( intArr, L1[1], R1[1], L2[1], R2[1] )
end

function integralSum_unsafe( intArr::AbstractArray{T,1}, 
                             L1::Int, 
                             R1::Int,
                             L2::Int, 
                             R2::Int
                           )::T where {T<:AbstractFloat}
    @inbounds begin  
        return  intArr[ R2+1 ] - intArr[ L2 ] - 
                intArr[ R1+1 ] + intArr[ L1 ] 
    end
end

""" 2D """
function integralSum_unsafe( intArr::AbstractArray{T,2}, 
                             TL1::Dims{2}, 
                             BR1::Dims{2},
                             TL2::Dims{2}, 
                             BR2::Dims{2} 
                           )::T where {T<:AbstractFloat}
    @inbounds begin  
        return intArr[TL2[1],TL2[2]] - intArr[BR2[1]+1,TL2[2]] - intArr[TL2[1],BR2[2]+1] + intArr[BR2[1]+1,BR2[2]+1] -
               intArr[TL1[1],TL1[2]] + intArr[BR1[1]+1,TL1[2]] + intArr[TL1[1],BR1[2]+1] - intArr[BR1[1]+1,BR1[2]+1]     
    end
end

""" 3D """
function integralSum_unsafe( intArr::AbstractArray{T,3}, 
                             TLF1::Dims{3}, 
                             BRB1::Dims{3},
                             TLF2::Dims{3}, 
                             BRB2::Dims{3}
                           )::T where {T<:AbstractFloat}
    @inbounds begin 
        return  intArr[BRB2[1]+1,BRB2[2]+1,BRB2[3]+1] - intArr[TLF2[1],TLF2[2],TLF2[3]] - intArr[TLF2[1],BRB2[2]+1,BRB2[3]+1] - intArr[BRB2[1]+1,TLF2[2],BRB2[3]+1] - intArr[BRB2[1]+1,BRB2[2]+1,TLF2[3]] + intArr[BRB2[1]+1,TLF2[2],TLF2[3]] + intArr[TLF2[1],BRB2[2]+1,TLF2[3]] + intArr[TLF2[1],TLF2[2],BRB2[3]+1] -
                intArr[BRB1[1]+1,BRB1[2]+1,BRB1[3]+1] + intArr[TLF1[1],TLF1[2],TLF1[3]] + intArr[TLF1[1],BRB1[2]+1,BRB1[3]+1] + intArr[BRB1[1]+1,TLF1[2],BRB1[3]+1] + intArr[BRB1[1]+1,BRB1[2]+1,TLF1[3]] - intArr[BRB1[1]+1,TLF1[2],TLF1[3]] - intArr[TLF1[1],BRB1[2]+1,TLF1[3]] - intArr[TLF1[1],TLF1[2],BRB1[3]+1]
    end    
end

# A version that is aware of "N", the size of the ROI after bound-safe cropping,
# and combined the integralSum and N according to the "two-input" function "op".

""" 1D """
function integralSumN_unsafe( intArr::AbstractArray{T,1}, 
                              L1::Int, 
                              R1::Int, 
                              L2::Int, 
                              R2::Int,
                              op::Function=(sum::T,n::T)->(sum/n)
                            )::T where {T<:AbstractFloat}
    @inbounds begin  
        return op( integralSum_unsafe( intArr, L1, R1, L2, R2 ), T( R2 - L2 + 1 - R1 + L1 - 1 ) )
    end
end

#
function integralSumN_unsafe( intArr::AbstractArray{T,1}, 
                              L1::Dims{1}, 
                              R1::Dims{1}, 
                              L2::Dims{1}, 
                              R2::Dims{1},
                              op::Function=(sum::T,n::T)->(sum/n)
                            )::T where {T<:AbstractFloat}
    return integralSumN_unsafe( intArr, L1[1], R1[1], L2[1], R2[1], op )
end

""" 2D """
function integralSumN_unsafe( intArr::AbstractArray{T,2}, 
                              TL1::Dims{2}, 
                              BR1::Dims{2}, 
                              TL2::Dims{2}, 
                              BR2::Dims{2},
                              op::Function=(sum::T,n::T)->(sum/n)
                            )::T where {T<:AbstractFloat}
    @inbounds begin  
        return op( integralSum_unsafe( intArr, TL1, BR1, TL2, BR2 ), T( prod( BR2 .- TL2 .+ 1 ) - prod( BR1 .- TL1 .+ 1) ) )
    end
end

""" 3D """
function integralSumN_unsafe( intArr::AbstractArray{T,3}, 
                              TLF1::Dims{3}, 
                              BRB1::Dims{3}, 
                              TLF2::Dims{3}, 
                              BRB2::Dims{3},
                              op::Function=(sum::T,n::T)->(sum/n)
                            )::T where {T<:AbstractFloat}
    @inbounds begin 
        return op( integralSum_unsafe( intArr, TLF1, BRB1, TLF2, BRB2, op ), T( prod( BRB2 .- TLF2 .+ 1 ) - prod( BRB1 .- TLF1 .+ 1 ) ) )
    end    
end