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

function localL2avg!( intA::IntegralArray{T,N},
                      intA2::IntegralArray{T,N},
                      output::AbstractArray{T,N},
                      rad::Dims{N} ) where {T<:AbstractFloat,N}

    localL2avg!( intA.arr, intA2.arr, output, rad )
    return nothing
end

"""
    the sum of all L2s between all pairs of pixels in a rectangular
    ROI is given by:

    sum_ij( (p_i - p_j)^2 ) = sum_ij( p_i^2 + p_j^2 - 2p_ip_j )
                            = Sij( p_i^2 ) + Sij( p_j^2 ) - Sij( 2p_ip_j )
                            = NjSi( p_i^2 ) + NiSj( p_j^2 ) - 2 Sij( p_ip_j )

    Nj = Ni = number of pixels in the ROI. Si( p_i^2 ) == Sj( p_j^2 ) are the
    sum of squared values within the ROI. Sij( p_ip_j ) is the sum of pairs of
    products between all pixels, which can be computed with Si( p_i ) * Sj( p_j ),
    where Si(p_i) == Sj(p_j) are the sum of pixels in the ROI.

    Thus, when computing the "sum of all pairs of L2s" within a ROI, most of the
    intermediate sums are the same and they can be combined into:
      
        2*S(ROI.^2) - 2S(ROI)^2

    When computing the "sum of all pairs of L2s" between two ROIs, the intermediate
    sums will be different and one has to compute the extended forumula:

        N2*S(ROI1.^2) + N1*S(ROI2.^2) - 2 * S(ROI1) * S(ROI2)

    The average difference is obtained by dividing by the product of the number
    of elements in both ROIS.
"""
function localL2avg!( intA::AbstractArray{T,N},
                      intA2::AbstractArray{T,N},
                      output::AbstractArray{T,N},
                      rad::Dims{N} ) where {T<:AbstractFloat,N}

    # 1-. output = sum( IA )
    localsums!( intA, output, rad ); 

    # 2-. output = -2 .* sum( IA )² 
    output .*= output;
    output .*= T(-2); 

    # 3-. output = 2 .* N .* sum( IA² ) - 2 .* sum( IA )²
    # TODO: remove allocations... partially removed by converting N to Float. 
    localsums_N_f_op!( intA2, output, rad, *, T(2), + )

    # 4-. output = ( 2 .* N .* sum( IA² ) - 2 .* sum( IA )² )/( N² )
    localN_op!( output, rad, (x)->(x*x), / )

    return nothing
end

# TODO: implementation with two radii
function localL2avg!( intA::AbstractArray{T,N},
                      intA2::AbstractArray{T,N},
                      output::AbstractArray{T,N},
                      tmp::AbstractArray{T,N},
                      rad1::Dims{N},
                      rad2::Dims{N} ) where {T<:AbstractFloat,N}

    @assert all( rad1 .<= rad2 )

    # output = sum( out )
    localsums!( intA, output, rad2 );

    # output = sum( out ) - sum( in ) = sum( ring )
    localsums_N_f_op!( intA, output, rad1, f=T(1), op=(out::T,res::T)->(out-res), Nop=(sum::T,n::T)->(sum) )

    # output = -2*sum( ring )*sum( in )
    localsums_N_f_op!( intA, output, rad1, f=T(-2), op=(out::T,res::T)->(out*res), Nop=(sum::T,n::T)->(sum) )

    # output = Nring*sum(in^2) - 2*sum(ring)*sum(ring)
    localsums!( intA2, tmp, rad1 )
    localN_op!( tmp, rad1, rad2, op=(out::T,n::T)->(out*n) )
    output .+= tmp; 

    # output = Nring*sum(in^2) + Nin*sum(ring^2) - 2*sum(ring)*sum(ring)
    localsums!( intA2, tmp, rad1, rad2 )
    localN_op!( tmp, rad1, op=(out::T,n::T)->(out*n) )
    output .+= tmp; 

    # 4-. output = ( Nring*sum(in^2) + Nin*sum(ring^2) - 2*sum(ring)*sum(ring) )/( Nin*Nring )
    localN_op!( output, rad1, op=(out::T,n::T)->(out/n) )
    localN_op!( output, rad1, rad2, op=(out::T,n::T)->(out/n) )

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
