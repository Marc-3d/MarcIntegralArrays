using BenchmarkTools
using MarcIntegralArrays


###### RECTANGULAR ROI

function integralSum_inputs( T::Type, dims::Vararg{Int} )
    TLF = rand.( UnitRange.( 1, dims .- 1 ) );
    BRB = rand.( UnitRange.( TLF .+ 1, dims  ) )
    return T, dims..., rand( T, dims .+ 1 ), TLF, BRB, T(1)
end

#### SAFE 
println(); 

# 1D
T, N, IA, L, R, f = integralSum_inputs( Float64, 100 ); 
safe_1D = @belapsed MarcIntegralArrays.integralSum( $IA, $L, $R, $f )

# 2D
T, Ny, Nx, IA, TL, BR, f = integralSum_inputs( Float64, 100, 100 ); 
safe_2D = @belapsed MarcIntegralArrays.integralSum( $IA, $TL, $BR, $f )

# 3D
T, Ny, Nx, Nz, IA, TLF, BRB, f = integralSum_inputs( Float64, 30, 30, 30 ); 
safe_3D = @belapsed MarcIntegralArrays.integralSum( $IA, $TLF, $BRB, $f )

println( "integralSum 1D ran in: $(safe_1D * 1000000000) ns" )
println( "integralSum 2D ran in: $(safe_2D * 1000000000) ns" )
println( "integralSum 3D ran in: $(safe_3D * 1000000000) ns" )

#### UNSAFE
println(); 

# 1D
T, N, IA, L, R, f = integralSum_inputs( Float64, 100 ); 
safe_1D = @belapsed MarcIntegralArrays.integralSum_unsafe( $IA, $L, $R, $f )

# 2D
T, Ny, Nx, IA, TL, BR, f = integralSum_inputs( Float64, 100, 100 ); 
safe_2D = @belapsed MarcIntegralArrays.integralSum_unsafe( $IA, $TL, $BR, $f )

# 3D
T, Ny, Nx, Nz, IA, TLF, BRB, f = integralSum_inputs( Float64, 30, 30, 30 ); 
safe_3D = @belapsed MarcIntegralArrays.integralSum_unsafe( $IA, $TLF, $BRB, $f )

println( "integralSum_unsafe 1D ran in: $(safe_1D * 1000000000) ns" )
println( "integralSum_unsafe 2D ran in: $(safe_2D * 1000000000) ns" )
println( "integralSum_unsafe 3D ran in: $(safe_3D * 1000000000) ns" )


###### RING-LIKE ROI

function integralSum_ring_inputs( T::Type, dims::Vararg{Int} )
    @assert all( dims .> 4 )
    # outer ring ( at least 4 units long on each dimension )
    TLF_out = rand.( UnitRange.( 1, dims .- 3 ) );
    BRB_out = rand.( UnitRange.( TLF_out .+ 3, dims  ) );
    # inner ring within outer ring
    TLF_in  = rand.( UnitRange.( TLF_out .+ 1, BRB_out .- 2 ) ); 
    BRB_in  = rand.( UnitRange.( TLF_in  .+ 1, BRB_out .- 1 ) );

    return T, dims..., rand( T, dims .+ 1 ), TLF_in, BRB_in, TLF_out, BRB_out, T(1)
end

#### SAFE 
println(); 

# 1D
T, N, IA, L_in, R_in, L_out, R_out, f = integralSum_ring_inputs( Float64, 100 ); 
safe_ring_1D = @belapsed MarcIntegralArrays.integralSum( $IA, $L_in, $R_in, $L_out, $R_out, $f )

# 2D
T, Ny, Nx, IA, TL_in, BR_in, TL_out, BR_out, f = integralSum_ring_inputs( Float64, 100, 100 ); 
safe_ring_2D = @belapsed MarcIntegralArrays.integralSum( $IA, $TL_in, $BR_in, $TL_out, $BR_out, $f )

# 3D
T, Ny, Nx, Nz, IA, TLF_in, BRB_in, TLF_out, BRB_out, f = integralSum_ring_inputs( Float64, 30, 30, 30 ); 
safe_ring_3D = @belapsed MarcIntegralArrays.integralSum( $IA, $TLF_in, $BRB_in, $TLF_out, $BRB_out, $f )

println( "integralSum 1D ran in: $(safe_ring_1D * 1000000000) ns" )
println( "integralSum 2D ran in: $(safe_ring_2D * 1000000000) ns" )
println( "integralSum 3D ran in: $(safe_ring_3D * 1000000000) ns" )
