function integralArraySQ( input::AbstractArray{T,N} 
                        ) where {T<:AbstractFloat,N}

    intArr = zeros( T, size(input) .+ 1 ); 
    integralArraySQ_unsafe!( intArr, input ); 
    return intArr
end
    
function integralArraySQ!( intArr::AbstractArray{T,N}, 
                           input::AbstractArray{T,N}
                          ) where {T,N}

    @assert size(intArr) == ( size(input) .+ 1 ); 
    integralArraySQ_unsafe!( intArr, input ); 
    return nothing    
end
		
function integralArraySQ_unsafe!( intArr::AbstractArray{T,2},
                                     img::AbstractArray{T,2}
                                 ) where {T<:Real}

    @inbounds for c in 2:size(img,2)+1
        cache = T(0.0)
        for r in 2:size(img,1)+1
           cache = img[r-1,c-1]*img[r-1,c-1] + cache + intArr[r,c-1] - intArr[r-1,c-1]
           intArr[r,c] = cache
        end
    end
    return nothing
end

    function integralArraySQ_unsafe!( intArr::AbstractArray{T,3},
		                         vol::AbstractArray{T,3}
		                     ) where {T<:Real}

        @inbounds for z in 2:size(vol,3)+1, c in 2:size(vol,2)+1
            cache = T(0.0)
            for r in 2:size(vol,1)+1
                cache = vol[r-1,c-1,z-1]*vol[r-1,c-1,z-1] + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1] + cache; 
                intArr[r,c,z] = cache
            end
        end
        return nothing
    end
