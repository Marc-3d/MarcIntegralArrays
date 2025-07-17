#### Singe-channel colors: Grayscale 

""" 
Array{Gray} -> IntegralArray

Computes integral array of grayscale input.
"""
function IntegralArray( 
    inp::AbstractArray{C,N};
    T = Float64,
    fun = (x)->(x)
) where {
    C<:ColorTypes.Color{<:Any,1},
    N
}
    IA = IntegralArray( T, size(inp) .+ 1 )
    integralArray!( IA, inp, fun=fun )

    return IA
end

"""
!IntegralArray, Array{Gray} -> nothing

Computes in-place integral array of a grayscale input
"""
function integralArray!( 
    intArr::IntegralArray{T,N}, 
    input::AbstractArray{C,N}; 
    fun::Function=(v)->(v)
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,1},
    N
}
     @assert size(intArr.arr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr.arr, input, fun ); 
     return nothing    
end

"""
!IntegralVector, Vector{Gray} -> nothing

1D single-channel in-place integral array computation
"""
function integralArray_unsafe!( 
    intArr::AbstractArray{T,1},
    vec::AbstractArray{C,1},
    fun::Function=(v)->(v)
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,1}
}
    cache = T(0)
    @inbounds for r in 2:size(vec,1)+1
          cache += fun(T(vec[r-1].val))
          intArr[r] = cache
    end
    return nothing
end

"""
!IntegralImage, Image{Gray} -> nothing

2D single-channel in-place integral array computation
"""
function integralArray_unsafe!( 
    intArr::AbstractArray{T,2},
    img::AbstractArray{C,2},
    fun::Function=(v)->(v) 
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,1}
}
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache += fun(T(img[r-1,c-1].val)) + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

# TODO: Implement 3D integral array

#### Multi-channel Colors: RGB, HSV, etc... 

"""
Array{Color{<:Any,3}}, channel -> IntegralArray

Computes integral array of multi-channel input from the desired channel
"""
function IntegralArray( 
    inp::AbstractArray{C,N};
    T::Type = Float64,
    channel::Int = 1
) where {
    C<:ColorTypes.Color{<:Any,3},
    N
}
    IA = IntegralArray( T, size(inp) .+ 1 )
    integralArray!( IA, inp, channel=channel )
    return IA
end

"""
IntegralArray (output), Array{Color{3}}, channel_id -> nothing

Computes in-place integral array of a multi-channel input from the desired channel
"""
function integralArray!( 
    intArr::IntegralArray{T,N}, 
    input::AbstractArray{C,N}; 
    fun::Function=(v)->(v),
    channel::Int = 1
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,3},
    N
}
     @assert size(intArr.arr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr.arr, input, fun, channel ); 
     return nothing    
end

"""
1D multi-channel in-place integral array computation
"""
function integralArray_unsafe!( 
    intArr::AbstractArray{T,1},
    vec::AbstractArray{C,1},
    fun::Function=(v)->(v),
    channel::Int = 1
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,3}
}
    cache = T(0)
    @inbounds for r in 2:size(vec,1)+1
          cache += fun( T( getfield( vec[r-1], channel ) ) )
          intArr[r] = cache
    end
    return nothing
end

"""
2D multi-channel in-place integral array computation
"""
function integralArray_unsafe!( 
    intArr::AbstractArray{T,2},
    img::AbstractArray{C,2},
    fun::Function=(v)->(v),
    channel::Int=1
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,3}
}
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache += fun( T( getfield( img[r-1,c-1], channel ) ) ) + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

# TODO: implement 3D case
