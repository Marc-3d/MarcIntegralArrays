# This is a very specific class for computing "integral dot products", which I first applied for 
# segmenting "contraction waves" in PIV vectorfields of the development of Tribolium Castaneum.
# TODO: properly define the integral dot product.


#=
    NOTE: Usually NC == ND. For instance, 2D vector fields are expressed by two (NC=2) two-dimensional 
    (ND=2) matrices: U & V. Similarly, 3D vector fields are expressed by three (NC=3) three-dimensional
    (ND=3) matrices: U, V & W. 
    
    This changes if we are dealing with datasets with time, such as 2D+t and 3D+t microscopy recordings. 
    In these cases, one could have:
    ( NC = 2, ND = 3 ): two 3D arrays, containing the spatiotemporal distribution of 2D vectors.
    ( NC = 3, ND = 4 ): three 4D arrays, containing the spatiotemporal distribution of 3D vectors.
    ( NC = 2, ND = 1 ): two 1D arrays, containing the timeseries for a single 2D vector.
    ( NC = 3, ND = 1 ): three 1D arrays, containing the timeseries for a single 3D vector.
=#
mutable struct IntegralVectorField{T,NC,ND}
    components::NTuple{NC,IntegralArray{T,ND}}
    magnitudes::IntegralArray{T,ND}
end

"""
    Transforming a vector field into an integral vector field, where an integral array is
    constructed for each component, in addition to an integral array of the magnitudes of
    the vector field.

    intVF2D = IntegralVectorField( U, V )
    intVF3D = IntegralVectorField( U, V, W )
"""
function IntegralVectorField( 
    components::Vararg{AbstractArray{R,N}}; 
    T = Float64
) where {
    R <: Real,
    N
}
    M = compute_magnitude( components... )
    integralVF = IntegralVectorField( IntegralArray.(components, T), IntegralArray(M, T) )
    return integralVF
end

Base.size( IVF::IntegralVectorField{T,NC,ND} ) where {T,NC,ND} = size( IVF.magnitudes.arr ) .- 1

# NOTE: lengths below only defined for 1D vector fields, aka timeseries of 2D or 3D vectors
""" length of the underlying data, aka length of the integral array - 1 """
Base.length( IVF::IntegralVectorField{T,NC,1} ) where {T,NC} = length( IVF.magnitudes.arr ) - 1 
""" length of the integral arrays, aka length of the data + 1 """
ilength( IVF::IntegralVectorField{T,NC,1} )  where {T,NC} = length( IVF.magnitudes.arr )

### UTILS

"""
    Average vectors

    input < IntegralVectorField( U, V[, W ] )
    output > Array( #Components, size_vectorfield... )
"""
function averageVF( 
    IVF::IntegralVectorField{T,NC,ND},
    rad::NTuple{ND,Int} 
) where {
    T,
    NC,
    ND
}
    avgVF = zeros( T, NC, size(IVF)... )
    for c in CartesianIndices( size(IVF) )
        TLF,BRB = c.I .- rad, c.I .+ rad
        for i in 1:NC
            avgVF[i,c.I...] = integralAverage( IVF.components[i].arr, TLF, BRB )
        end
    end
    return avgVF
end

"""
    Collectiveness: 
        indot of vector at "p" with vectors around it.
"""
function collectiveness( 
    IVF::IntegralVectorField{T,NC,ND};
    rad::NTuple{ND,Int} 
) where {
    T,
    NC,
    ND
}
    out = zeros( T, size(IVF) )
    for c in CartesianIndices( out )
        TLF,BRB = c.I .- rad, c.I .+ rad
        ndot = indot( IVF, UnitRange.( c.I, c.I ), UnitRange.( TLF, BRB ) ); 
        score = T(0.0)
        # for i in 1:NC
        #     s1 = IVF.components[i][ c.I, c.I ] 
        #     s2 = IVF.components[i][ TLF, BRB ] - s1
        #     score += s1 * s2
        # end
        # m1 = IVF.magnitudes[ c.I, c.I ]
        # m2 = IVF.magnitudes[ TLF, BRB ] - m1
        out[ c ] = ndot # ( score )/( m1 * m2 )
    end
    return out
end


"""
    The function below computes the "integral normalized dot product" (indot, for short) between
    two ROIs of an `IntegralVectorField`. Each ROI can be given as:
        1D: "one UnitRange" ...or... "a tuple of two ints" 
        ND: "a 2-tuple of NTuples{ND,Int} ...or... or an NTuple{ND,UnitRange}

    1D: indot( IVF, (1,5), (4,8) ) -- ROI1 = 1:5, ROI2 = 4:8
    1D: indot( IVF, 1:5, 4:8 ) -- ROI1 = 1:5, ROI2 = 4:8
    ND: indot( IVF, ((2,3,4),(8,10,6)), ((3,5,1),(4,5,3)) ) -- ROI1 = (2:8,3:10,4:6), ROI2=(3:4,5:5,1:3)
    ND: indot( IVF, (2:8, 3:10, 4:6), (3:4, 5:5, 1:3) ) -- ROI1 = (2:8,3:10,4:6), ROI2=(3:4,5:5,1:3)
"""
function indot( 
    IVF::IntegralVectorField{T,NC,ND}, 
    ROI1::Union{NTuple{2,NTuple{ND,Int}},NTuple{ND,UnitRange}},
    ROI2::Union{NTuple{2,NTuple{ND,Int}},NTuple{ND,UnitRange}} 
) where {
    T,
    NC,
    ND
}
    indot = T(0.0)
    for i in 1:NC
        indot += IVF.components[i][ ROI1... ] * IVF.components[i][ ROI2... ]
    end
    indot /= IVF.magnitudes[ ROI1... ] * IVF.magnitudes[ ROI2... ]
    return indot
end

function indots( 
    IVF::IntegralVectorField{T,NC,ND}, 
    radii::Dims{ND}=Tuple(zeros(Int,ND))
) where {
    T,
    NC,
    ND
}
    indots = zeros( T, size( IVF ) )
    for c in CartesianIndices( size( IVF ) )
        indots[c] = indot( IVF, (Tuple(c).-radii,Tuple(c).+radii), (Tuple(c).-radii,Tuple(c).+radii)  )
    end
    return indots
end

function indots_( 
    IVF::IntegralVectorField{T,NC,ND};
    left_1::Dims{ND}=Tuple(zeros(Int,ND)),
    right_1::Dims{ND}=Tuple(zeros(Int,ND)),
    off_1::Dims{ND}=Tuple(zeros(Int,ND)),
    left_2::Dims{ND}=Tuple(zeros(Int,ND)),
    right_2::Dims{ND}=Tuple(zeros(Int,ND)),
    off_2::Dims{ND}=Tuple(zeros(Int,ND))
) where {
    T,
    NC,
    ND
}
    indots = zeros( T, size( IVF ) )
    indots!( 
        indots,
        IVF,
        left_1=left_1,
        right_1=right_1,
        off_1=off_1,
        left_2=left_2,
        right_2=right_2,
        off_2=off_2
    )
    return indots
end

function indots!( 
    out::AbstractArray{T,ND},
    IVF::IntegralVectorField{T,NC,ND};
    left_1::Dims{ND}=Tuple(zeros(Int,ND)),
    right_1::Dims{ND}=Tuple(zeros(Int,ND)),
    off_1::Dims{ND}=Tuple(zeros(Int,ND)),
    left_2::Dims{ND}=Tuple(zeros(Int,ND)),
    right_2::Dims{ND}=Tuple(zeros(Int,ND)),
    off_2::Dims{ND}=Tuple(zeros(Int,ND))
) where {
    T,
    NC,
    ND
}
    for c in CartesianIndices( size( IVF ) )
        out[c] = indot( 
            IVF, 
            (Tuple(c).+off_1.+left_1,Tuple(c).+off_1.+right_1), 
            (Tuple(c).+off_2.+left_2,Tuple(c).+off_2.+right_2)  
        )
    end
    return nothing
end

###

function indot( 
    IVF::IntegralVectorField{T,NC,1}, 
    ROI1::NTuple{2,Int},
    ROI2::NTuple{2,Int} 
) where {
    T,
    NC
}
    indot = T(0)
    for i in 1:NC
        indot += IVF.components[i][ ROI1... ] * IVF.components[i][ ROI2... ]
    end
    indot /= IVF.magnitudes[ ROI1... ] * IVF.magnitudes[ ROI2... ]
    return indot
end

function indot( 
    IVF::IntegralVectorField{T,NC,1}, 
    ROI1::UnitRange,
    ROI2::UnitRange 
) where {
    T,
    NC
}
    indot = T(0)
    for i in 1:NC
        indot += IVF.components[i][ ROI1 ] * IVF.components[i][ ROI2 ]
    end
    indot /= IVF.magnitudes[ ROI1 ] * IVF.magnitudes[ ROI2 ]
    return indot
end


"""
    The function below computes the "integral normalized dot product" (indot, for short) between
    one ROI of an `IntegralVectorField` and a reference vector. 
"""
function indot( 
    IVF::IntegralVectorField{T,NC,ND}, 
    ROI1::Union{NTuple{2,NTuple{ND,Int}},NTuple{ND,UnitRange}},
    refV::NTuple{NC,T},
    magV::T 
) where {
    T,
    NC,
    ND
}
    indot = T(0.0)
    for i in 1:NC
        indot += IVF.components[i][ ROI1... ] * refV[i]
    end
    indot /= IVF.magnitudes[ ROI1... ] * magV
    return indot
end

function indots( 
    IVF::IntegralVectorField{T,NC,ND}, 
    refV::NTuple{NC,T},
    radii::Dims{ND}=Tuple(zeros(Int,ND))
) where {
    T,
    NC,
    ND
}
    magV = sqrt( sum( refV .* refV ) )
    indots = zeros( T, size( IVF ) )
    for c in CartesianIndices( size( IVF ) )
        indots[c] = indot( IVF, (Tuple(c).-radii,Tuple(c).+radii), refV, magV  )
    end
    return indots
end

function indots_( 
    IVF::IntegralVectorField{T,NC,ND}, 
    refV::NTuple{NC,T}, 
    magV=sqrt(sum(refV.*refV));
    left::Dims{ND}=Tuple(zeros(Int,ND)),
    right::Dims{ND}=Tuple(zeros(Int,ND))
) where {
    T,
    NC,
    ND
}
    indots = zeros( T, size( IVF ) )
    indots!( 
        indots, 
        IVF, 
        refV, 
        magV,
        left=left, 
        right=right
    )
    return indots
end


function indots!( 
    out::AbstractArray{T,ND},
    IVF::IntegralVectorField{T,NC,ND}, 
    refV::NTuple{NC,T}, 
    magV=sqrt(sum(refV .*refV) );
    left::Dims{ND}=Tuple(zeros(Int,ND)),
    right::Dims{ND}=Tuple(zeros(Int,ND))
) where {
    T,
    NC,
    ND
}
    for c in CartesianIndices( size( IVF ) )
        out[c] = indot( IVF, (Tuple(c).+left,Tuple(c).+right), refV, magV  )
    end
    return nothing
end

###

function inmags( 
    IVF::IntegralVectorField{T,NC,ND}, 
    ROI1::Union{NTuple{2,NTuple{ND,Int}},NTuple{ND,UnitRange}},
) where {
    T,
    NC,
    ND
}
    return IVF.magnitudes[ ROI1... ]
end

function inmags( 
    IVF::IntegralVectorField{T,NC,ND};
    left::Dims{ND}=Tuple(zeros(Int,ND)),
    right::Dims{ND}=Tuple(zeros(Int,ND))
) where {
    T,
    NC,
    ND
}
    indots = zeros( T, size( IVF ) )
    for c in CartesianIndices( size( IVF ) )
        indots[c] = inmags( IVF, (Tuple(c).+left,Tuple(c).+right)  )
    end
    return indots
end


####################################

function sliding_idnot( IVF::IntegralVectorField{T,NC,1}, 
                        Nprev::Int ) where {T,NC}

    N   = length( IVF ); 
    out = zeros( T, N ); 
    for i in 1:N
        out[i] = indot( IVF, i:i, i-Nprev:i-1 ) # i think this is bound-safe...
    end
    return out
end

####

"""
    this class keeps N entries, and stores the "parent" of each entry. The way that
    it is implemented, if two entries are combined, their parents will be set to the
    smallest of both entries. The specific functions of DisjointMinSets are defined
    further below in the file. 
"""
struct DisjointMinSets
    parents::Vector{Int}
    DisjointMinSets(n::Integer) = new([1:n;]) # [1:n;] == collect( 1:n )
end
DisjointMinSets() = DisjointMinSets(0)

"""
    This function computes the sliding_idnot at a desired scale, and merges vectors in
    a "connected-component-style" if the idnot score is higher than a certain threshold.

    In the initial application of `_label_components!`, Albl is expected to be a vector
    of labels where each element has a different label, which can be its own linear index.
    Similarly, sets should initially be set to: sets = DisjointMinSets( length(Albl) ). 

    In further applications of `_label_components!`, Albl and sets will be reused. 
"""
function _label_components!( Albl::AbstractArray{Int,1},
                             IVF::IntegralVectorField{T,NC,1},
                             Nprev::Int = 3,
                             dot_th::T = T(cosd(20)),
                             sets::DisjointMinSets = DisjointMinSets(length(Albl))
                           ) where {NC,T<:AbstractFloat}

    @assert length( Albl ) == length( IVF )

    score = sliding_idnot( IVF, Nprev ); 

    for i in 2:length(score)

        label = Albl[i]

        if score[i] > dot_th 
            for ii in max(1,i-Nprev):i
                newlabel = Albl[ ii ]
                if label != newlabel
                    label = union!(sets, label, newlabel)  # ...merge labels...
                else
                    label = newlabel  # ...and assign the smaller to current pixel
                end
            end
        end

        Albl[i] = label
    end

    # Now parse sets to find the labels
    newlabel = minlabel(sets)
    for i = 1:length(Albl)
        Albl[i] = newlabel[find_root!(sets, Albl[i])]
    end

    return Albl
end

"""
    Out-of-place implementation, it accepts multiple scales like so:

    labels = _label_components( IVF, scale1, scale2, scale3, dot_th=cosd(15) )
"""
function _label_components( IVF::IntegralVectorField{T,NC,1},
                            Nprevs::Int...;
                            dot_th::T = T(cosd(20)),
                          ) where {NC,T<:AbstractFloat}
    N = length( IVF );
    Albl = collect( 1:N );
    sets = DisjointMinSets( N );
    for n in Nprevs
        _label_components!( Albl, IVF, n, dot_th, sets )
    end
    return Albl
end

function _label_components_mask( IVF::IntegralVectorField{T,NC,1},
                                 Nprevs::Int...;
                                 dot_th::T = T(cosd(20)),
                               ) where {NC,T<:AbstractFloat}

    lbls = _label_components( IVF, Nprevs..., dot_th=dot_th );
    mask = zeros( UInt8, length( lbls ) )
    val  = 0;  
    for i in 2:length( lbls )
        val = ( lbls[i] == lbls[i-1] ) ? val : 1 - val
        mask[i] = val
    end
    return reshape( mask, (1,length(mask)) );
end

#####

function compute_magnitude( 
    components::Vararg{AbstractArray{T,N}} 
) where {
    T,
    N
}
    M = zeros( T, size(components[1]) ); 
    @inbounds for i in 1:length(components)
        M .+= components[i] .^ 2 
    end
    @inbounds for i in 1:length(M)
        M[i] = sqrt( M[i] )
    end
    return M 
end

# Functions copied from github.com/Marc-3d/directionCCL, which is in turn an adaptation of
# of ImageComponentAnalysis.jl. Suprisingly, these 50 lines of code are 50% of what is 
# needed to implement connected components labelling. 

import Base.push!  

function find_root!(sets::DisjointMinSets, m::Integer)
    p = sets.parents[m]   # don't use @inbounds here, it might not be safe
    @inbounds if sets.parents[p] != p
        sets.parents[m] = p = find_root_unsafe!(sets, p)
    end
    p
end

# an unsafe variant of the above
function find_root_unsafe!(sets::DisjointMinSets, m::Int)
    @inbounds p = sets.parents[m]
    @inbounds if sets.parents[p] != p
        sets.parents[m] = p = find_root_unsafe!(sets, p)
    end
    p
end

function union!(sets::DisjointMinSets, m::Integer, n::Integer)
    mp = find_root!(sets, m)
    np = find_root!(sets, n)
    if mp < np
        sets.parents[np] = mp
        return mp
    elseif np < mp
        sets.parents[mp] = np
        return np
    end
    mp
end

function push!(sets::DisjointMinSets)
    m = length(sets.parents) + 1
    push!(sets.parents, m)
    m
end

function minlabel(sets::DisjointMinSets)
    out = Vector{Int}(undef, length(sets.parents))
    k = 0
    for i = 1:length(sets.parents)
        if sets.parents[i] == i
            k += 1
        end
        out[i] = k
    end
    out
end

