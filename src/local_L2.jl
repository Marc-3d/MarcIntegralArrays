"""
    For each position in the input image, this function computes "the sum of all pairs of L2 differences between all pixels within a rectangular ROI of half-size 'rad'". For each position, this quantity is computed as: 

    \$ ( 2 .* N .* sum( I² ) - 2 .* sum( I )² )/( N² ) 

    In other words, from "the sum of squared pixel values within the square ROI (time 2N)" we substract "(2 times) the sum of pixel within the square ROI, squared". This quantity can be computed very efficiently with two integral arrays and in-place integral sums.    

"""
function localL2avg( 
    img::AbstractArray{T,N}, 
    rad ::Dims{N}
) where {
    T,
    N
}
   intAL2 = IntegralArraysL2( img )
   out = localL2avg( intAL2, rad )
   return out  
end

function localL2avg( 
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    output = zeros( T, size(intAL2.IA.arr) .- 1 ); 
    localL2avg!( output, intAL2, rad ) 
    return output
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2.IA, intAL2.IA2, rad )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    rad::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intA.arr, intA2.arr, rad )
    return nothing
end

#=
    The sum of all L2s between all pairs of pixels in a rectangular ROI is given by:

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
function localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    rad::Dims{N};
    op::Function=(out::T,ret::T)->(ret*ret)::T
) where {
    T<:AbstractFloat,
    N
}
    # 1-. output = -2 .* sum( IA )²
    localsums!( 
        output, 
        intA, 
        rad,
        op=op
    ); 
    output .*= T(-2)

    # 2-. output = 2 .* N .* sum( IA² ) - 2 .* sum( IA )²
    localsumNs!( 
        output, 
        intA2, 
        rad, 
        f=T(2), 
        op=+
    )

    # 3-. output = ( 2 .* N .* sum( IA² ) - 2 .* sum( IA )² )/( N² )
    localN!( 
        output, 
        rad, 
        op=(o::T,n::T)->(o/(n*n)) 
    )

    return nothing
end

"""
    This function computes the "sum of all pairs of L2 differences" (like the function above), but it does so in a ring-like ROI, instead of a ROI. Basically, a ring-like ROI can be computed by considering a the sum of value in a large square... and substracting the sum of values within a smaller square around the center of the larger square.
"""

function localL2avg( 
    input::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    intAL2 = IntegralArraysL2( input ); 
    output = zeros( eltype( intAL2.IA.arr ),size( input ) ); 
    localL2avg!( output, intAL2, rad_in, rad_out )
    return output
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2::IntegralArraysL2{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2.IA, intAL2.IA2, rad_in, rad_out )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intA.arr, intA2.arr, rad_in, rad_out )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N}; 
    op::Function=(out::T,in::T)->(in*in)::T
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

    # 3-. output = ( 2 .* N .* sum( ring² ) - 2 .* sum( ring )² )/( N² )
    localN!( 
        output, 
        rad_in,
        rad_out, 
        op=(o::T,n::T)->(o/(n*n)) 
    )

    return nothing
end

######

function localL2avg( 
    input::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    intAL2 = IntegralArraysL2( input ); 
    output = zeros( eltype( intAL2.IA.arr ),size( input ) ); 
    tmp    = zeros( eltype( intAL2.IA.arr ),size( input ) );
    localL2avg!( output, intAL2, tmp, rad_in, rad_mid, rad_out )
    return output 
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2::IntegralArraysL2{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2.IA, intAL2.IA2, tmp, rad_in, rad_mid, rad_out )
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
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intA.arr, intA2.arr, tmp, rad_in, rad_mid, rad_out )
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
    op::Function=(out::T,in::T)->(in)::T
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

    # 4-. output = ( Nring*sum(in^2) + Nin*sum(ring^2) - 2*sum(ring)*sum(ring) )/( Nin*Nring )
    localN!( 
        output, 
        rad_in, 
        op=(out::T,n::T)->(out/n)
    )
    localN_op!( 
        output, 
        rad_mid, 
        rad_out, 
        op=(out::T,n::T)->(out/n)
    )

    return nothing
end


################################                       

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

