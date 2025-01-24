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

#####

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

"""
    The sum of all L2s between all pairs of pixels in a rectangular ROI is given by:

    sum_ij( (p_i - p_j)^2 ) 
        = sum_ij( p_i^2 + p_j^2 - 2p_ip_j )
        = Sij( p_i^2 ) + Sij( p_j^2 ) - Sij( 2p_ip_j )
        = NjSi( p_i^2 ) + NiSj( p_j^2 ) - 2 Sij( p_ip_j )

    Nj == Ni = number of pixels in the ROI. Si( p_i^2 ) == Sj( p_j^2 ) are the sum of squared values within the ROI. Sij( p_ip_j ) is the sum of pairs of products between all pixels, which can be computed with Si( p_i ) * Sj( p_j ), where Si( p_i ) == Sj( p_j ) are the sum of pixels in the ROI.

    Thus, when computing the "sum of all pairs of L2s" within a ROI, most of the intermediate sums are the same and they can be combined into:
      
        2*S(ROI.^2) - 2S(ROI)^2

    When computing the "sum of all pairs of L2s" between two ROIs, the intermediate sums will be different and one has to compute the extended forumula:

        N2*S(ROI1.^2) + N1*S(ROI2.^2) - 2 * S(ROI1) * S(ROI2)

    The average difference is obtained by dividing by the product of the number of elements in both ROIS.
"""
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

# TODO: implementation with two radii
"""

"""

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2::IntegralArraysL2{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2.IA, intAL2.IA2, tmp, rad_in, rad_out, op=op )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intA.arr, intA2.arr, tmp, rad_in, rad_out, op=op )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N}; 
    op::Function=(out::T,in::T)->(out+in)::T
) where {
    T<:AbstractFloat,
    N
}
    @assert all( rad_in .<= rad_out )

    # output = sum( out ) - sum( in ) = sum( ring )
    localsums!( 
        output, 
        intA,
        rad_in,
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
    # ...

    # output = Nring*sum(in^2) - 2*sum(ring)*sum(ring)
    localsums!( 
        tmp,
        intA2, 
        rad_in,
        op=op
    )
    localN!( 
        tmp, 
        rad_in, 
        rad_out, 
        op=(out::T,n::T)->(out*n)::T
    )
    output .+= tmp; 

    # output = Nring*sum(in^2) + Nin*sum(ring^2) - 2*sum(ring)*sum(ring)
    localsums!( 
        tmp, 
        intA2, 
        rad_in, 
        rad_out,
        op=(o::T,r::T)->(r)::T
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
        rad_in, 
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

