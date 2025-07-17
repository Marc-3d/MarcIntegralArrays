#### THESE FUNCTIONS ARE LIKE THE ONES ABOVE; BUT THEY ALSO TAKE A MASK
# TODO: test all masked functions          

"""
    square ROI WITH ITSELF
"""
function localL2avg( 
    img::AbstractArray{C,N}, 
    mask::AbstractArray{I,N},
    rad ::Dims{N}=Tuple(ones(Int,N)).*3;
    T::Type=Float64
) where {
    C<:Union{Real,Color{<:Any,1}},
    I<:Real,
    N
}
   intAL2M = IntegralArraysL2M( img, mask )
   out = localL2avg( intAL2M, rad )
   return out  
end

function localL2avg( 
    intAL2M::IntegralArraysL2M{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    output = zeros( T, size(intAL2M.IA.arr) .- 1 ); 
    localL2avg!( output, intAL2M, rad ) 
    return output
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2M::IntegralArraysL2M{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2M.IA, intAL2M.IA2, intAL2M.IAM, rad )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    intAM::IntegralArray{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}
    _localL2avg!( output, intA.arr, intA2.arr, intAM.arr, rad )
    return nothing
end

#=
    We assume that intA and intA2 have been computed on a maske input, where everything outside the mask are zeros.

    We merely need to replace the "N"s by the actual sum of mask elements in each ROI.
=#
function _localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    intAM::AbstractArray{T,N},
    rad::Dims{N};
    op::Function=(out::T,ret::T)->(out-2*ret*ret)::T
) where {
    T<:AbstractFloat,
    N
}

    # 1-. output = 2 .* sum( IA² )
    localsums!( 
        output, 
        intA2, 
        rad, 
        f=T(2), 
        op=(out::T,ret::T)->(ret)::T
    )

    # 2-. output = 2 .* N .* sum( IA² )
    localsums!( 
        output, 
        intAM, 
        rad, 
        op=(out::T,ret::T)->(out*ret)::T
    )

    # 3-. output = 2 .* N .* sum( IA² ) - 2 .* sum( IA )²
    localsums!( 
        output, 
        intA, 
        rad,
        op=(out::T,ret::T)->(out-2*ret*ret)::T
    ); 

    # 4-. output = ( 2 .* N .* sum( IA² ) - 2 .* sum( IA )² )/( N² )
    localsums!( 
        output, 
        intAM,
        rad, 
        op=(out::T,ret::T)->(out/(ret*ret))::T
    )

    return nothing
end

"""

"""

function localL2inter( 
    img::AbstractArray{C,N}, 
    mask::AbstractArray{I,N},
    rad ::Dims{N}=Tuple(ones(Int,N)).*3;
    T::Type=Float64
) where {
    C<:Union{Real,Color{<:Any,1}},
    I<:Real,
    N
}
   intAL2M_M = IntegralArraysL2M( img .* mask, mask )
   intAL2M_m = IntegralArraysL2M( img .* ( mask .< 0.5 ), mask .< 0.5 )
   out = zeros( T, size( img ) ); 
   tmp = zeros( T, size( img ) ); 

   _localL2inter!( out, tmp, intAL2M_M, intAL2M_m, rad )

   return out  
end


# intra differences: masked (M) vs unmasked (m)
function _localL2inter!( 
    output::AbstractArray{T,N},
    tmp::AbstractArray{T,N},
    intAL2M_M::IntegralArraysL2M{T,N},
    intAL2M_m::IntegralArraysL2M{T,N},
    rad::Dims{N}
) where {
    T<:AbstractFloat,
    N
}

    # 1.1-. output = sum( ( IA .* mask )² )
    localsums!( 
        output, 
        intAL2M_M.IA2.arr, 
        rad, 
        f=T(1), 
        op=(out::T,ret::T)->(ret)::T
    )

    # 1.2-. output = N_m * sum( ( IA .* mask )² )
    localsums!( 
        output, 
        intAL2M_m.IAM.arr, 
        rad, 
        f=T(1), 
        op=(out::T,ret::T)-> (out*ret)::T
    )

    # 2.1-. tmp = sum( ( IA .* !mask )² )
    localsums!( 
        tmp, 
        intAL2M_m.IA2.arr, 
        rad, 
        f=T(1), 
        op=(out::T,ret::T)->(ret)::T
    )

    # 2.2-. tmp = N_M * sum( ( IA .* !mask )² )
    localsums!( 
        tmp, 
        intAL2M_M.IAM.arr, 
        rad, 
        f=T(1), 
        op=(out::T,ret::T)-> (out * ret)::T
    )

    # 3-. output = N_m * sum( ( IA .* mask )² ) + N_M * sum( ( IA .* !mask )² )
    output .+= tmp; 

    # 4.1-. tmp = sum( ( IA .* mask ) )
    localsums!( 
        tmp, 
        intAL2M_M.IA.arr, 
        rad, 
        op=(out::T,ret::T)->(ret)::T
    )

    # 4.2-. tmp = sum( ( IA .* mask ) ) * sum( ( IA .* !mask ) )
    localsums!( 
        tmp, 
        intAL2M_m.IA.arr, 
        rad, 
        op=(out::T,ret::T)->(out * ret)::T
    )

    # 5-. output = N_m * sum( ( IA .* !mask )² ) + N_M * sum( ( IA .* !mask )² ) + sum( ( IA .* mask ) ) * sum( ( IA .* !mask ) )
    output .+= tmp; 

    # 6.1-. output = ( N_m * sum( ( IA .* !mask )² ) + N_M * sum( ( IA .* !mask )² ) + sum( ( IA .* mask ) ) * sum( ( IA .* !mask ) ) ) / ( N_m )
    localsums!( 
        output, 
        intAL2M_m.IAM.arr, 
        rad,
        op=(out::T,ret::T)->(out/ret)::T
    ); 


    # 6.2-. output = ( N_m * sum( ( IA .* !mask )² ) + N_M * sum( ( IA .* !mask )² ) + sum( ( IA .* mask ) ) * sum( ( IA .* !mask ) ) ) / ( N_m * N_M )
    localsums!( 
        output, 
        intAL2M_M.IAM.arr, 
        rad,
        op=(out::T,ret::T)->(out/ret)::T
    ); 

    return nothing
end





"""
    ring-like ROI masked L2
"""

function localL2avg( 
    input::AbstractArray{T,N},
    mask::AbstractArray{I,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    I,
    N
}
    intAL2M = IntegralArraysL2M( input, mask ); 
    output  = zeros( eltype( intAL2M.IA.arr ),size( input ) ); 
    localL2avg!( output, intAL2M, rad_in, rad_out )
    return output
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2M::IntegralArraysL2M{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2M.IA, intAL2M.IA2, intAL2M.IAM, rad_in, rad_out )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    intAM::IntegralArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intA.arr, intA2.arr, intAM.arr, rad_in, rad_out )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    intAM::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_out::Dims{N}; 
    op::Function=(out::T,in::T)->(in*in)::T
) where {
    T<:AbstractFloat,
    N
}
    @assert all( rad_in .<= rad_out )

    # 2-. output = 2 .* sum( ring² )
    localsums!( 
        output, 
        intA2, 
        rad_in,
        rad_out,
        f=T(2), 
        op=(out::T,ret::T)->(ret)::T
    )

    # 2-. output = 2 .* N .* sum( ring² )
    localsums!( 
        output, 
        intAM, 
        rad_in,
        rad_out,
        op=(out::T,ret::T)->(out * ret)::T
    )

    # 3-. output = 2 .* N .* sum( ring² ) - 2 * sum( ring )^2
    localsums!( 
        output, 
        intA,
        rad_in,
        rad_out,
        op=(out::T,ret::T)->(out-2*ret*ret)::T
    );

    # 4-. output = ( 2 .* N .* sum( ring² ) - 2 .* sum( ring )² )/( N² )
    localsums!( 
        output, 
        intM,
        rad_in,
        rad_out,
        op=(out::T,ret::T)->(out/(ret*ret))::T
    );

    return nothing
end

######

function localL2avg( 
    input::AbstractArray{T,N},
    mask::AbstractArray{I,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    I,
    N
}
    intAL2M = IntegralArraysL2M( input, mask ); 
    output  = zeros( eltype( intAL2M.IA.arr ),size( input ) ); 
    tmp     = zeros( eltype( intAL2M.IA.arr ),size( input ) );
    localL2avg!( output, intAL2M, tmp, rad_in, rad_mid, rad_out )
    return output 
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intAL2M::IntegralArraysL2M{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intAL2M.IA, intAL2M.IA2, intAL2M.IAM, tmp, rad_in, rad_mid, rad_out )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::IntegralArray{T,N},
    intA2::IntegralArray{T,N},
    intAM::IntegralArray{T,N},
    tmp::AbstractArray{T,N},
    rad_in::Dims{N},
    rad_mid::Dims{N},
    rad_out::Dims{N};
) where {
    T<:AbstractFloat,
    N
}
    localL2avg!( output, intA.arr, intA2.arr, intAM.arr, tmp, rad_in, rad_mid, rad_out )
    return nothing
end

function localL2avg!( 
    output::AbstractArray{T,N},
    intA::AbstractArray{T,N},
    intA2::AbstractArray{T,N},
    intAM::AbstractArray{T,N},
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
        op=(out::T,ret::T)->(ret)::T
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
        op=(out::T,ret::T)->(ret)::T
    )
    localsums!( 
        tmp,
        intAM, 
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
        op=(out::T,ret::T)->(ret)::T
    )
    localsums!( 
        tmp,
        intAM, 
        rad_in,
        op=*
    )
    output .+= tmp; 

    # 4-. output = ( Nring*sum(in^2) + Nin*sum(ring^2) - 2*sum(ring)*sum(ring) )/( Nin*Nring )
    localsums!( 
        tmp,
        intAM, 
        rad_in,
        op=(out::T,n::T)->(out/n)::T
    )
    localsums!( 
        tmp,
        intAM, 
        rad_mid,
        rad_out,
        op=(out::T,n::T)->(out/n)::T
    )

    return nothing
end