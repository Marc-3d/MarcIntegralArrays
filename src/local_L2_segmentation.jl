# Taken from my experiments with tribolium castaneum

"""
    This is class that holds two integral arrays for any input array.
    Namely, it hold a "standard integral array", and a "square
    integral array". Both of these are required to compute local L2s.
"""
mutable struct IntegralArraysL2{T,N}
    IA::IntegralArray{T,N}
    IA2::IntegralArray{T,N}
end

function integralArraysL2( inp::AbstractArray{T,N} ) where {T,N}

    IA_size = size(inp) .+ 1 
    out = IntegralArraysL2{T,N}( IntegralArray( zeros( T, IA_size ) ), 
                                 IntegralArray( zeros( T, IA_size ) ) )

    integralArray!( out.IA, inp )
    integralArray!( out.IA2, inp, (x)::T->(x*x) )

    return out
end

function integralArraysL2!( IAL2::IntegralArraysL2{T,N}, inp::AbstractArray{T,N} ) where {T,N}

    integralArray!( IAL2.IA, inp )
    integralArray!( IAL2.IA2, inp, (x)::T->(x*x) )
    return nothing
end

####################

function sqr!( arr::AbstractArray )
    @inbounds @fastmath arr .= sqrt.( arr )
    return nothing
end

"""
    For each pixel, we compute local (dis)similarity measures:
        pixel - local mean 
        avgL2 around pixel
        std of avgL2 
"""
function local_L2_stat_segmentation( inp::AbstractArray{T,N};
                                     avg_r::Dims{N}=Tuple(ones(Int,N)),
                                     dif_r::Dims{N}=Tuple(ones(Int,N)),
                                     std_r::Dims{N}=Tuple(ones(Int,N)),
                                     avg_f=1,
                                     std_f=1 ) where {T,N}

    IAL2 = integralArraysL2( inp  ); 
    tmp1 = zeros(  T  , size(inp) ); 
    tmp2 = zeros(  T  , size(inp) ); 
    mask =  ones( Bool, size(inp) ); 

    # tmp1 = inp
    tmp1 .= inp; 
    # tmp1 = ( inp - avg ) = ( inp - sum(I)/N_ )
    localsums_N_f_op!( IAL2.IA.arr, tmp1, avg_r, f=T(-1) )

    # tmp2 = avg(L2(I))    
    localL2avg!( IAL2.IA, IAL2.IA2, tmp2, dif_r )
    sqr!( tmp2 )

    ###
    localL2avg!( IAL2.IA, IAL2.IA2, tmp1, avg_r )

    IAL2.IA.arr  .= 0.0
    integralArray!( IAL2.IA, tmp2 );
    localsums_N_f_op!( IAL2.IA.arr, tmp2, std_r )

    return tmp1 ./ tmp2
    ###
    
    # mask1 = local average dif > global average dif.
    avg_dif = T( sum( tmp1 ) * avg_f / length( tmp1 ) ); 
    mask .= tmp1 .> tmp2; 
    #
    IAL2.IA.arr  .= 0.0
    IAL2.IA2.arr .= 0.0
    integralArray!( IAL2.IA , tmp2, );
    integralArray!( IAL2.IA2, tmp2, (x)::T->(x*x)); 
    localL2avg!( IAL2.IA, IAL2.IA2, tmp2, std_r )
    sqr!( tmp2 )

    avg_std = T( sum( tmp1 ) * avg_f / length( tmp1 ) ); 

    mask .*= tmp1 .> ( tmp2 .* std_f )
    mask .*= tmp2 .> avg_std

    return tmp2
end

function local_L2_stat_segmentation( inp::AbstractArray{T,N};
                                     avg_r::Dims{N}=Tuple(ones(Int,N)),
                                     dif_r::Dims{N}=Tuple(ones(Int,N)),
                                     std_r::Dims{N}=Tuple(ones(Int,N)),
                                     avg_f=1,
                                     std_f=1 ) where {T,N}

    IAL2 = integralArraysL2( inp  ); 
    tmp1 = zeros(  T  , size(inp) ); 
    tmp2 = zeros(  T  , size(inp) ); 
    mask =  ones( Bool, size(inp) ); 

    # Sum of differences of pixel - ROI
    localL2avg!( IAL2.IA, IAL2.IA2, tmp1, avg_r )
    sqr!( tmp1 )

    # tmp2 = avg(L2(I))    
    localL2avg!( IAL2.IA, IAL2.IA2, tmp2, dif_r )
    sqr!( tmp2 )

    IAL2.IA.arr  .= 0.0
    integralArray!( IAL2.IA, tmp2 );
    localsums_N_f_op!( IAL2.IA.arr, tmp2, std_r )

    return tmp1 ./ tmp2
end

function local_L2_tmp_data( inp::AbstractArray{T,N} ) where {T,N}
    IAL2 = integralArraysL2( inp  ); 
    tmp1 = zeros(  T  , size(inp) ); 
    tmp2 = zeros(  T  , size(inp) ); 
    mask =  ones( Bool, size(inp) );
    return IAL2, tmp1, tmp2, mask
end

function local_L2_stat_segmentation_2!( inp::AbstractArray{T,N},
                                        IAL2::IntegralArraysL2{T,N},
                                        tmp1::AbstractArray{T,N},
                                        tmp2::AbstractArray{T,N},
                                        mask::AbstractArray{Bool,N};
                                        rad0::Dims{N}=Tuple( ones(Int,N)),
                                        rad1::Dims{N}=Tuple(zeros(Int,N)),
                                        rad2::Dims{N}=Tuple( ones(Int,N)),
                                        rad3::Dims{N}=Tuple( ones(Int,N)),
                                        rad4::Dims{N}=Tuple( ones(Int,N)),
                                        avg_f=1,
                                        std_f=1 ) where {T,N}

    # Local & global intensity based thresholding

    localsums_N_f_op!( IAL2.IA.arr, tmp1, rad0, op=(out::T,res::T)->(res), Nop=(sum::T,n::T)->(sum/n) )
    global_th = sum( tmp1 )/length( tmp1 )

    mask .*= inp .> tmp1
    mask .*= inp .> global_th

    # Segmenting regions where the L2 difs are a larger scale (GS) are greater than the avg L2 dif at low scales (LS).
    # This is a sign that there are structures of thickness (LS), so that locally the average L2 difs at LS are low...
    # and the average L2 difs at the larger scale will be high.

    # tmp1 = sum of differences of pixel - ROI
    tmp1 .= 0.0
    localL2avg!( IAL2.IA.arr, IAL2.IA2.arr, tmp1, tmp2, rad1, rad2 )
    sqr!( tmp1 )

    # tmp2 = avg(L2(I))    
    tmp2 .= 0.0
    localL2avg!( IAL2.IA, IAL2.IA2, tmp2, rad3 )
    sqr!( tmp2 )

    # tmp2 = local( avg(L2(I)) )
    IAL2.IA.arr .= 0.0
    integralArray!( IAL2.IA, tmp2 );
    localsums_N_f_op!( IAL2.IA.arr, tmp2, rad4, op=(out::T,res::T)->(res), Nop=(sum::T,n::T)->(sum/n) )

    mask .*= ( tmp1 ./ tmp2 ) .> avg_f

    return nothing
end
