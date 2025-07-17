################################                       

"""   
    Computes the standard deviation in a local square ROI around each pixel. 
"""
function localSTD( 
    img::AbstractArray{C,N}, 
    rad ::Dims{N}; 
    T=Float64
) where {
    C<:Union{Real,Color{<:Any,1}},
    N
}
   intAL2 = IntegralArraysL2( img, T )
   out = localSTD( intAL2, rad )
   return out  
end

#####

function localSTD( 
    intAL2::IntegralArraysL2{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    output = zeros( T, size(intAL2.IA.arr) .- 1 ); 
    localSTD!( output, intAL2, rad ) 
    return output
end

"""
    The sum of all L2s between all pixels and the mean intensity in a rectangular ROI is given by:

    var_{ROI] = sum_j( (I_j - M)^2 ) / N_j
        = sum_j( I_j^2 + M^2 - 2MI_j ) / N_j
        = ( sum_j( I_j^2 ) + sum_j( M^2 ) - sum_j( 2MI_j ) ) / N_j 
        = ( sum_j( I_j^2 ) + N_J * M^2 - 2 * M * sum_j( I_j ) ) / N_j 
        = ( sum_j( I_j^2 ) + N_j * sum_j( I_j )^2/N_j^2 - 2 * M * sum_j( I_j ) ) / N_j
        = ( sum_j( I_j^2 ) + sum_j( I_j )^2/N_j - 2 * sum_j( I_j )^2/N_j ) / N_j
        = ( sum_j( I_j^2 ) - sum_j( I_j )^2/N_j ) / N_j
"""
function localSTD!( 
    output::AbstractArray{T,N},
    intL2::IntegralArraysL2{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    # 1-. output = - sum( IA )^2 / N
    localsums!( 
        output, 
        intL2.IA.arr, 
        rad,
        op=(o::T,r::T)->(r*r) 
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
        intL2.IA2.arr, 
        rad, 
        op=(o::T,r::T) -> (o+r)
    )

    # 3-. output = ( sum( IA² ) - sum( IA )^2 / N )/( N )
    localN!( 
        output, 
        rad, 
        op=(o::T,n::T)->(o/n) 
    )

    return nothing
end