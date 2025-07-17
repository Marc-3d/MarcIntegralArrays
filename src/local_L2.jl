"""
    This function computes "the sum of all pairwise L2 differences between all pixels within a rectangular ROI", or "L2avg". This operation is applied locally around each position in the input array. Namely, the function considers a rectangular ROI of half-side 'rad' around each position: UnitRange( pos .- rad, pos .+ rad ). From each ROI, the L2avg is computed from: 

    ``( 2 .* N .* \\sum_R( I[ x ]² ) - 2 .* \\sum_R( I[ x ] )² )/( N² )``

    In other words, from "the sum of squared pixel values within the square ROI (time 2N)" we substract "(2 times) the sum of pixel within the square ROI, squared". This quantity can be computed very efficiently with two integral arrays and in-place integral sums.    
"""

# Real numbers and grayscale images
function localL2avg( 
    img::AbstractArray{C,N}, 
    rad::Dims{N}=Tuple(ones(Int,N)).*3; 
    T::Type=Float64,
    average=true
) where {
    C<:Union{Real,Color{<:Any,1}},
    N
}
    intAL2 = IntegralArraysL2( img, T )
    out = localL2avg( intAL2, rad, average=average )
    return out  
end

# Multi-channel images: RGB, HSV, ... 
function localL2avg( 
    img::AbstractArray{C,N}, 
    rad::Dims{N}=Tuple(ones(Int,N)).*3; 
    T::Type=Float64, 
    average=true
) where {
    C<:Color{<:Any,3},
    N
}
    intAL2 = IntegralArraysL2( img, T, field=1, average )
    out = localL2avg( intAL2, rad, average=average )
    tmp = zeros( T, size( out ) ) 

    for c in 2:3
       tmp .= 0.0
       integralArraysL2!( intAL2, img, c )
       localL2avg!( tmp, intAL2, rad, average=average )
       out .+= tmp
    end
    return out  
end

#-----------------------------

function localL2avg( 
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}; 
    average=true
) where {
    T<:AbstractFloat,
    N
}
    output = zeros( T, size(intAL2.IA.arr) .- 1 ); 
    localL2avg!( output, intAL2, rad, average=average ) 
    return output
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}; 
    average=true
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2.IA, intAL2.IA2, rad, average=average )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    rad::Dims{N};
    average=true
) where {
    T<:AbstractFloat,
    N
}
    _localL2avg!( output, intA.arr, intA2.arr, rad, average=average )
    return nothing
end

#=
    The sum of all pairwise L2 differences between all pairs of pixels in a rectangular ROI is given by:

    \$ sum_i{ sum_j{ ( p_i - p_j )^2 } }
    = sum_i{ sum_j{ ( p_i^2 + p_j^2 - 2p_ip_j ) } }
    = sum_i{ sum_j{ p_i^2 } } + sum_i{ sum_j{  p_j^2 } } - sum_i{ sum_j{  2p_ip_j } }
    = N_j * sum_i{ p_i^2 } + N_i * sum_j{  p_j^2 } - 2 sum_i{ sum_j{  p_ip_j } }

    where, N_j == N_i are number of pixels in the ROI. 
    where, sum_i{ p_i^2 } and sum_j{ p_j^2 } are the sum of squared values within the ROI. 
    where, sum_i{ sum_j{  p_ip_j } } is the sum of pairs of products between all pixels, which can be computed with sum_i( p_i ) * sum_j( p_j ). 

    In our case, sum_i{ p_i^2 } == sum_j{ p_j^2 } and sum_i( p_i ) == sum_j( p_j ) since we are computing the sum of differences within the same ROI. Thus, when computing the "sum of all pairs of L2s" within a signel ROI, most of the intermediate sums are the same and they can be combined into:
      
    \$ 2 * N * sum_ROI{ p.^2 } - 2 sum_ROI{ p }^2

    When computing the "sum of all pairs of L2s" between two ROIs (R1 and R2), the intermediate sums will be different and one has to compute the extended forumula:

    \$ N_R2 * sum_R1{ p .^2 } + N_R1 * sum_R2{ p .^2 } - 2 * sum_R1{ p } * sum_R2{ p }

    The average difference is obtained by dividing by the product of the number of elements in both ROIS.
=#
function _localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    rad::Dims{N};
    op::Function=(out::T,ret::T)->(ret*ret)::T,
    average=true
) where {
    T<:AbstractFloat,
    N
}
    # 1-. output = -2 .* sum( IA )²
    localsums!( 
        output, 
        intA, 
        rad,
        op=(o::T,r::T)->(r*r) 
    ); 
    output .*= T(-2)

    # 2-. output = 2 .* N .* sum( IA² ) - 2 .* sum( IA )²
    localsumNs!( 
        output, 
        intA2, 
        rad, 
        f=T(2), 
        op=(o::T,r::T)->(o+r) 
    )   

    if !average
        return nothing
    end

    # 3-. output = ( 2 .* N .* sum( IA² ) - 2 .* sum( IA )² )/( N² )
    localN!( 
        output, 
        rad, 
        op=(o::T,n::T)->(o/(n*n)) 
    )

    return nothing
end

"""
    This function computes the "sum of all pairs of L2 differences" (like the function above), but it does so in a ring-like ROI. Basically, a ring-like ROI can be computed by considering a the sum of value in a large square... and substracting the sum of values within a smaller square around the center of the larger square.
"""

function localL2avg( 
    input::AbstractArray{C,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
    T::Type=Float64, 
    average=true
) where {
    C<:Union{Real,Color{<:Any,1}},
    N
}
    intAL2 = IntegralArraysL2( input, T ); 
    output = zeros( T,size( input ) ); 
    localL2avg!( output, intAL2, rad_in, rad_out, average=average )
    return output
end

function localL2avg( 
    input::AbstractArray{C,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
    T::Type=Float64,
    average=true
) where {
    C<:Color{<:Any,3},
    N
}
    intAL2 = IntegralArraysL2( input, T, field=1 )
    output = zeros( T, size( input ) )
    localL2avg!( output, intAL2, rad_in, rad_out, average=average )
    
    tmp = zeros( T, size( output ) ) 
    for c in 2:3
       tmp .= 0.0
       integralArraysL2!( intAL2, input, c )
       localL2avg!( tmp, intAL2, rad_in, rad_out, average=average )
       output .+= tmp
    end
    return output  
end

#-----------------------------

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2::IntegralArraysL2{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
    average=true
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2.IA, intAL2.IA2, rad_in, rad_out, average=average )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
    average=true
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intA.arr, intA2.arr, rad_in, rad_out, average=average )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N}; 
    op::Function=(out::T,in::T)->(in*in)::T, 
    average=true
) where {
    T<:AbstractFloat,
    N
}
    @assert all( rad_in .<= rad_out )

    # 1-. output = -2 * sum( ring )^2
    localsums!( 
        output, 
        intA,
        rad_in,
        rad_out,
        op=op
    );
    output .*= T(-2)

    # 2-. output = 2 .* N .* sum( ring² ) - 2 .* sum( ring )²
    localsumNs!( 
        output, 
        intA2, 
        rad_in,
        rad_out,
        f=T(2), 
        op=+
    )

    if !average
        return nothing
    end

    # 3-. output = ( 2 .* N .* sum( ring² ) - 2 .* sum( ring )² )/( N² )
    localN!( 
        output, 
        rad_in,
        rad_out, 
        op=(o::T,n::T)->(o/(n*n)) 
    )

    return nothing
end

"""
    This function computes the "sum of all pairs of L2 differences" (like the function above) between an inner square ROI, and an outer ring-like ROI. The inner squre ROI is defined by a single "rad" parameters, while the outer ring-like ROI is defined by two "rad" parameters.  
"""
function localL2avg( 
    input::AbstractArray{C,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
    T::Type=Float64, 
    average=true
) where {
    C<:Union{Real,Color{<:Any,1}},
    N
}
    intAL2 = IntegralArraysL2( input, T ); 
    output = zeros( T, size( input ) ); 
    tmp    = zeros( T, size( input ) );
    localL2avg!( output, intAL2, tmp, rad_in, rad_mid, rad_out, average=average )
    return output 
end

function localL2avg( 
    input::AbstractArray{C,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
    T::Type=Float64, 
    average=true
) where {
    C<:Color{<:Any,3},
    N
}
    intAL2 = IntegralArraysL2( input, T ); 
    output = zeros( T, size( input ) ); 
    tmp1   = zeros( T, size( input ) );
    tmp2   = zeros( T, size( input ) ); 
    localL2avg!( output, intAL2, tmp1, rad_in, rad_mid, rad_out, average=average )

    for c in 2:3
       tmp1 .= 0.0
       tmp2 .= 0.0
       integralArraysL2!( intAL2, input, c )
       localL2avg!( tmp1, intAL2, tmp2, rad_in, rad_mid, rad_out, average=average )
       output .+= tmp1
    end  

    return output 
end

 
#--------------------------------

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2::IntegralArraysL2{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
    average=true
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2.IA, intAL2.IA2, tmp, rad_in, rad_mid, rad_out, average=average )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
    average=true
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intA.arr, intA2.arr, tmp, rad_in, rad_mid, rad_out, average=average )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N}; 
    op::Function=(out::T,in::T)->(in)::T, 
    average=true
) where {
    T<:AbstractFloat,
    N
}
    @assert all( rad_in .<= rad_out )

    # output = sum( ring )
    localsums!( 
        output, 
        intA,
        rad_mid,
        rad_out,
        op=op
    );

    # output = -2*sum( ring )*sum( in )
    localsums!( 
        output, 
        intA,
        rad_in,
        f=T(-2),
        op=*
    )

    # output = Nring*sum(in^2) - 2*sum(ring)*sum(in)
    localsums!( 
        tmp,
        intA2, 
        rad_in,
        op=op
    )
    localN!( 
        tmp, 
        rad_mid, 
        rad_out, 
        op=*
    )
    output .+= tmp; 

    # output = Nring*sum(in^2) + Nin*sum(ring^2) - 2*sum(ring)*sum(ring)
    localsums!( 
        tmp, 
        intA2, 
        rad_mid, 
        rad_out,
        op=op
    )
    localN!( 
        tmp, 
        rad_in, 
        op=(out::T,n::T)->(out*n)
    )
    output .+= tmp; 

    if !average
        return nothing
    end

    # 4-. output = ( Nring*sum(in^2) + Nin*sum(ring^2) - 2*sum(ring)*sum(ring) )/( Nin*Nring )
    localN!( 
        output, 
        rad_in, 
        op=(out::T,n::T)->(out/n)
    )
    localN!( 
        output, 
        rad_mid, 
        rad_out, 
        op=(out::T,n::T)->(out/n)
    )

    return nothing
end



#### THESE FUNCTIONS COMPUTE THE L2 DIFFERENCES BETWEEN ONE SIGNLE PIXEL (NOT A ROI) AND ITS SURROUNDING PIXELS                   

"""   
"""
function localL2avg_( 
    img::AbstractArray{T,N}, 
    rad ::Dims{N}
) where {
    T,
    N
}
   intAL2 = IntegralArraysL2( img )
   out = localL2avg_( img, intAL2, rad )
   return out  
end

#####

function localL2avg_( 
    img::AbstractArray{T,N}, 
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    output = zeros( T, size(img) ); 
    intAL2 = IntegralArraysL2( img )
    localL2avg_!( output, img, intAL2, rad ) 
    return output
end

function localL2avg_( 
    img::AbstractArray{T,N}, 
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    output = zeros( T, size(img) ); 
    localL2avg_!( output, img, intAL2, rad ) 
    return output
end

function localL2avg_!( 
    output::AbstractArray{T,N},
    img::AbstractArray{T,N}, 
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    localL2avg_!( output, img, intAL2.IA, intAL2.IA2, rad )
    return nothing
end

function localL2avg_!( 
    output::AbstractArray{T,N},
    img::AbstractArray{T,N}, 
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    rad::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg_!( output, img, intA.arr, intA2.arr, rad )
    return nothing
end

"""
    The sum of all L2s between all pairs of pixels in a rectangular ROI is given by:

    sum_j( (I - p_j)^2 ) / Nj
        = sum_j( I^2 + p_j^2 - 2Ip_j ) / Nj
        = ( Sj( I^2 ) + Sj( p_j^2 ) - Sj( 2Ip_j ) ) / Nj 
        = ( Nj( I^2 ) + Sj( p_j^2 ) - 2ISj( p_j ) ) / Nj
        = I^2 + ( Sj( p_j^2 ) - 2ISj( p_j ) ) / Nj
"""
function localL2avg_!( 
    output::AbstractArray{T,N},
    img::AbstractArray{T,N}, 
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    rad::Dims{N};
    op::Function=(out::T,ret::T)->(ret)::T
) where {
    T<:AbstractFloat,
    N
}

    # 1-. output = -2 .* I .* sum( IA )
    localsums!( 
        output, 
        intA, 
        rad,
        f=T(-2),
        op=op
    ); 
    output .*= img

    # 2-. output = sum( IA² ) - 2 .* I .* sum( IA )
    localsums!( 
        output, 
        intA2, 
        rad, 
        op=+
    )

    # 3-. output = ( sum( IA² ) - 2 .* I .* sum( IA ) )/( N )
    localN!( 
        output, 
        rad, 
        op=(o::T,n::T)->(o/n) 
    )

    # 3-. output = I.^2 .+ ( sum( IA² ) - 2 .* I .* sum( IA )² )/( N )
    output .+= img .^ 2

    return nothing
end

################################                       

"""   
"""
function localSTD_( 
    img::AbstractArray{T,N}, 
    rad ::Dims{N}
) where {
    T,
    N
}
   intAL2 = IntegralArraysL2( img )
   out = localSTD_( img, intAL2, rad )
   return out  
end

#####

function localSTD_( 
    img::AbstractArray{T,N}, 
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    output = zeros( T, size(img) ); 
    localSTD_!( output, img, intAL2, rad ) 
    return output
end

function localSTD_!( 
    output::AbstractArray{T,N},
    img::AbstractArray{T,N}, 
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    localSTD_!( output, img, intAL2.IA, intAL2.IA2, rad )
    return nothing
end

function localSTD_!( 
    output::AbstractArray{T,N},
    img::AbstractArray{T,N}, 
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    rad::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localSTD_!( output, img, intA.arr, intA2.arr, rad )
    return nothing
end

"""
    The sum of all L2s between all pairs of pixels in a rectangular ROI is given by:

    sum_j( (I_j - M)^2 ) / Nj
        = sum_j( I_j^2 + M^2 - 2MI_j ) / Nj
        = ( Sj( I_j^2 ) + Sj( M^2 ) - Sj( 2MI_j ) ) / Nj 
        = ( Sj( I_j^2 ) + NjSj( I_j )^2/Nj^2 - 2MSj( I_j ) ) / Nj
        = ( Sj( I_j^2 ) + Sj( I_j )^2/Nj - 2Sj( I_j )^2/Nj ) / Nj
        = ( Sj( I_j^2 ) - Sj( I_j )^2/Nj ) / Nj
"""
function localSTD_!( 
    output::AbstractArray{T,N},
    img::AbstractArray{T,N}, 
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    rad::Dims{N};
    op::Function=(out::T,ret::T)->(ret*ret)::T
) where {
    T<:AbstractFloat,
    N
}

    # 1-. output = - sum( IA )^2 / N
    localsums!( 
        output, 
        intA, 
        rad,
        op=op
    ); 
    output .*= -1
    localN!( 
        output, 
        rad, 
        op=(o::T,n::T)->(o/n) 
    )

    # 2-. output = sum( IA² ) - sum( IA )^2 / N
    localsums!( 
        output, 
        intA2, 
        rad, 
        op=+
    )

    # 3-. output = ( sum( IA² ) - sum( IA )^2 / N )/( N )
    localN!( 
        output, 
        rad, 
        op=(o::T,n::T)->(o/n) 
    )

    return nothing
end
