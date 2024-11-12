# TODO: optimize this
"""
     For each position in the input, this function computes "the
    sum of all pairs of L2 differences between the pixels within
    a rectangular ROI of half-size 'rad'". For each position, this
    quantity is computed as: 

         ( 2 .* N .* sum( I² ) - 2 .* sum( I )² )/( N² ) 

    In other words, from "the local sum of squared values" (times 2N)
    we substract "the local sum of values, squared" (times 2) for each
    position. This is can be computed very efficiently with two 
    integral arrays and in-place integral sums.    
"""
function localL2avg( img, rad )
   
   intA  = integralArray( img ); 
   intA2 = integralArraySQ( img ); 
   
   sumIA = zeros( eltype(img), size(img) ); 
   localsums!( intA, sumIA, rad ); 

   sumIA2 = zeros( eltype(img), size(img) ); 
   localsums!( intA2, sumIA2, rad );

   N = localN( img, rad ); 

   return ( 2 .* N .* sumIA2 .- 2 .* sumIA .* sumIA )./( N.^2 )  
end

function localL2avg!( intA::AbstractArray{T,N},
                      intA2::AbstractArray{T,N},
                      output::AbstractArray{T,N},
                      rad::Dims{N} ) where {T<:AbstractFloat,N}

    # 1-. output = sum( IA )
    localsums!( intA, output, rad ); 

    # 2-. output = -2 .* sum( IA )² 
    output .*= output
    output .*= T(-2); 

    # 3-. output = 2 .* N .* sum( IA² ) - 2 .* sum( IA )²
    # TODO: remove allocations... partially removed by converting N to Float. 
    localsums_N_f_op!( intA2, output, rad, *, T(2), + )

    # 4-. output = ( 2 .* N .* sum( IA² ) - 2 .* sum( IA )² )/( N² )
    localN_op!( output, rad, (x)->(x*x), / )

    return nothing

end
                       

# TODO: properify
function fun_sketch( img;
                     dm=10, 
                     r1=(5,5), 
                     r2=(2,2),
                     r3=(2,2),
                     r4=(10,10),
                     fmax=2, ovp=(0,0), f=1 )

    E = local_L2avg( img, r1 ); 
    E .*= E .> dm

    mn, mx = local_extrema( E, r2, fmax=fmax, f=f ); 

    mF = Float64.( mx ); 
    vF = mF .* img;  
    mB = Float64.( mn ); 
    vB = mB .* img; 

    # dif FB > dif BB
    NF  = localsums( mF, r3 );
    NB  = localsums( mB, r4 );
    sF  = localsums( vF, r3 );
    sB  = localsums( vB, r4 );
    sF2 = localsums( vF.^2, r3 );
    sB2 = localsums( vB.^2, r4 );
    difFB = ( NB .* sF2 .+ NF .* sB2 .- 2 .* sF .* sB )./( NB.*NF )
    difBB = ( 2 .* NB .* sB2 .- 2 .* sB .* sB )./(NB.*NB)

    return mF .* ( difFB .> difBB ), E, mx, NF, NB, sF, sB, sF2, sB2, difFB, difBB
end
