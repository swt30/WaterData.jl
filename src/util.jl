# Utility functions used throughout the package

using Iterators: partition


# Custom exceptions

"Error raised when a particular method needs to be implemented for some type"
immutable NotImplementedError <: Exception; end


# Duplicate-finding

""" Get indices of the unique items in `itr`
    For non-unique items, returns the highest index.

    ```jl
    julia> WaterData.unique_indices([1,1,4,7,4,4])
    3-element Array{Int64,1}:
     2
     6
     4
    ```
    """
unique_indices(itr) = indexin(collect(unique(itr)), collect(itr))


# Miscellaneous utility funcs

typealias VectorPair{T} NTuple{2, Vector{T}}

""" Map a function `f` across rows of the array `A`

    ```jl
    julia> A = [1 2 3;
                4 5 6];
    julia> WaterData.maprows(sum, A)
    2x1 Array{Int64,2}:
      6
     15
    ```
    """
maprows(f, A::Matrix) = mapslices(f, A, 2)

""" Get an iterator of each adjacent pair in a sequence

    ```jl
    julia> seq = [1,2,3,4];
    julia> collect(WaterData.adjacentpairs(seq))
    3-element Array{Tuple{Int64,Int64},1}:
     (1,2)
     (2,3)
     (3,4)
    ```
    """
adjacentpairs(itr) = partition(itr, 2, 1)

"2D analogue of the cross product, *v* × *w* = v₁w₂ - v₂w₁"
crossprod2d(v::Vector, w::Vector) = v[1]*w[2] - v[2]*w[1]

"Is x outside the closed interval [a, b]?"
isoutside(x, a, b) = x < a || x > b

"Check if the line segments `AB` and `CD` cross"
function intersects(A::Vector, B::Vector, C::Vector, D::Vector)
    # AB is the line from p to p + r
    # CD is the line from q to q + s
    p = A
    r = B - A
    q = C
    s = D - C

    rxs = crossprod2d(r, s)
    qpxr = crossprod2d((q - p), r)
    t = crossprod2d((q - p), s) ./ rxs
    u = qpxr ./ rxs

    if (rxs ≠ 0) && (0 <= t <= 1) && (0 <= u <= 1)
        # cross or touch
        return true
    else
        # may be colinear or disjoint, but they don't cross
        return false
    end
end
intersects(AB::VectorPair, CD::VectorPair) = intersects(AB[1], AB[2], CD[1], CD[2])
