# Grayscale images (1D Color)

function IntegralArray( 
    inp::AbstractArray{C,N},
    T = Float64
) where {
    C<:ColorTypes.Color{<:Any,1},
    N
}
    IA = IntegralArray( T, size(inp) .+ 1 )
    integralArray!( IA, inp, (x)->(T(x)) )

    return IA
end

function integralArray!( 
    intArr::IntegralArray{T,N}, 
    input::AbstractArray{C,N}, 
    fun::Function=(v)->(T(v))
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,1},
    N
}
     @assert size(intArr.arr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr.arr, input, fun ); 
     return nothing    
end

function integralArray_unsafe!( 
    intArr::AbstractArray{T,1},
    vec::AbstractArray{C,1},
    fun::Function=(v)->(T(v))
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,1}
}
    cache = T(0)
    @inbounds for r in 2:size(vec,1)+1
          cache += fun(vec[r-1].val)
          intArr[r] = cache
    end
    return nothing
end

function integralArray_unsafe!( 
    intArr::AbstractArray{T,2},
    img::AbstractArray{C,2},
    fun::Function=(v)->(T(v)) 
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,1}
}
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache += fun(img[r-1,c-1].val) + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

# Multi-channel Colors: RGB, HSV, etc... 

function IntegralArray( 
    inp::AbstractArray{C,N},
    T = Float64;
    field::Int=1
) where {
    C<:ColorTypes.Color{<:Any,3},
    N
}
    IA = IntegralArray( T, size(inp) .+ 1 )
    integralArray!( IA, inp, (x)->(T(x)), field )

    return IA
end

function integralArray!( 
    intArr::IntegralArray{T,N}, 
    input::AbstractArray{C,N}, 
    fun::Function=(v)->(T(v)),
    field::Int = 1
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,3},
    N
}
     @assert size(intArr.arr) == ( size(input) .+ 1 ); 
     integralArray_unsafe!( intArr.arr, input, fun, field ); 
     return nothing    
end

function integralArray_unsafe!( 
    intArr::AbstractArray{T,1},
    vec::AbstractArray{C,1},
    fun::Function=(v)->(T(v)),
    field::Int = 1
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,3}
}
    cache = T(0)
    @inbounds for r in 2:size(vec,1)+1
          cache += fun( getfield( vec[r-1], field ) )
          intArr[r] = cache
    end
    return nothing
end

function integralArray_unsafe!( 
    intArr::AbstractArray{T,2},
    img::AbstractArray{C,2},
    fun::Function=(v)->(T(v)),
    field::Int=1
) where {
    T<:AbstractFloat,
    C<:ColorTypes.Color{<:Any,3}
}
    @inbounds for c in 2:size(img,2)+1
        cache = T(0)
        for r in 2:size(img,1)+1
            cache += fun( getfield( img[r-1,c-1], field ) ) + intArr[r,c-1] - intArr[r-1,c-1]
            intArr[r,c] = cache
        end
    end
    return nothing
end

