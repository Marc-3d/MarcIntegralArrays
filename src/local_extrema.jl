# TODO: import and adapt 3D version 

function local_extrema(
    img::Array{T,2},
    rad::Dims{2}; 
    avg_th=0,
    fmin=1, 
    fmax=1, 
    ovp=(0,0), 
    f=1 
) where {
    T<:Real
}
    minima = zeros(Bool,size(img)); 
    maxima = zeros(Bool,size(img)); 
    intA   = IntegralArray(img); 
    local_extrema!( 
        minima, 
        maxima,
        intA, 
        rad,
        avg_th=avg_th, 
        fmin=fmin, 
        fmax=fmax, 
        ovp=ovp, 
        f=f 
    );
    return minima, maxima
end

function local_extrema!( 
    minima::Array{Bool,2}, 
    maxima::Array{Bool,2}, 
    intA::IntegralArray{T,2},
    rad::Dims{2};
    avg_th=0,
    fmin=1, 
    fmax=1, 
    ovp=0, 
    f=1 
) where {
    T<:AbstractFloat
}    
    # convenient quantities
    isize = size(intA.arr) .- 1; 
    ystep = ( 2*rad[1] + 1 - ovp[1], 0  ); 
    xstep = ( 0, 2*rad[2] + 1  - ovp[2] );
    
    # these will accomodate the means and signs of the 9x9 mean neighbourhood
    sums9 = zeros(   T  , 3, 3 ); 
    signs = zeros( Int64, 3, 3 ); 

    for cidx in CartesianIndices( isize ) #  in 1:isize[2], y in 1:isize[1]
        
        TL = max.( cidx.I .- rad,   1   ); 
        BR = min.( cidx.I .+ rad, isize );
        NN = prod( BR .- TL .+ 1 ); 
        
        sum_ = integralSum_unsafe( intA, TL, BR );
        avg  = sum_/NN    
        ( avg <= avg_th ) && ( continue; )
        
        idx = 1;
        @inbounds for xoff in -1:1, yoff in -1:1
            off = ystep .* yoff .+ xstep .* xoff; 
            sums9[idx] = integralSum( intA, TL .+ off, BR .+ off ); 
            signs[idx] = sign( sum_*f - sums9[idx] ); 
            idx += 1;   
        end

        score_min = 0
        score_max = 0
        @inbounds for idx in 1:4
            score_min += ( signs[idx] == -1 && signs[end-idx+1] == -1 )   
            score_max += ( signs[idx] ==  1 && signs[end-idx+1] ==  1 )
        end
        
        @inbounds minima[ cidx ] = score_min >= fmin; 
        @inbounds maxima[ cidx ] = score_max >= fmax; 
    end

    return nothing
end
