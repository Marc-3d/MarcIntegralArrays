using MarcIntegralArrays, BenchmarkTools

#### Rectangular ROI
println()

# 1D
N = 100; T = Float32; IA = rand( T, N+1 ); output = zeros( T, N ); rad = (5,); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
localSums_1D = @belapsed MarcIntegralArrays.localsums!( $IA, $output, $rad, f=$f, op=$op )

println( "localSums 1D ran in: $(localSums_1D * 1000000000 / N) ns per element" ); 

# 2D
Ny, Nx = 100, 100; T = Float32; IA = rand( T, Ny+1, Nx+1 ); output = zeros( T, Ny, Nx ); rad = (5,5); f = T(1); op = (o,i)->(i);
localSums_2D = @belapsed MarcIntegralArrays.localsums!( $IA, $output, $rad, f=$f, op=$op )

println( "localSums 2D ran in: $(localSums_2D * 1000000000 / ( Ny * Nx ) ) ns per element" ); 

# 3D
Ny, Nx, Nz = 20, 20, 20; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); output = zeros( T, Ny, Nx, Nz ); rad = (5,5,5); f = T(1); op = (o,i)->(i);
localSums_3D = @belapsed MarcIntegralArrays.localsums!( $IA, $output, $rad, f=$f, op=$op )

println( "localSums 3D ran in: $(localSums_3D * 1000000000 / ( Ny * Nx * Nz ) ) ns per element" ); 


#### Ring ROI
println()

# 1D
N = 100; T = Float32; IA = rand( T, N+1 ); output = zeros( T, N ); rad_in = (5,); rad_out = (10,); f = Float32(1); op = (o::Float32,i::Float32)->(i)::Float32;
localSums_1D_ring = @belapsed MarcIntegralArrays.localsums!( $IA, $output, $rad_in, $rad_out, f=$f, op=$op )

println( "localSums 1D ran in: $(localSums_1D_ring * 1000000000 / N) ns per element" ); 

# 2D
Ny, Nx = 100, 100; T = Float32; IA = rand( T, Ny+1, Nx+1 ); output = zeros( T, Ny, Nx ); rad_in = (5,5); rad_out = (10,10);  f = T(1); op = (o,i)->(i);
localSums_2D_ring = @belapsed MarcIntegralArrays.localsums!( $IA, $output, $rad_in, $rad_out, f=$f, op=$op )

println( "localSums 2D ran in: $(localSums_2D_ring * 1000000000 / ( Ny * Nx ) ) ns per element" ); 

# 3D
Ny, Nx, Nz = 20, 20, 20; T = Float32; IA = rand( T, Ny+1, Nx+1, Nz+1 ); output = zeros( T, Ny, Nx, Nz ); rad_in = (5,5,5); rad_out = (10,10,10); f = T(1); op = (o,i)->(i);
localSums_3D_ring = @belapsed MarcIntegralArrays.localsums!( $IA, $output, $rad_in, $rad_out, f=$f, op=$op )

println( "localSums 3D ran in: $(localSums_3D_ring * 1000000000 / ( Ny * Nx * Nz ) ) ns per element" ); 
