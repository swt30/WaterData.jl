using VoronoiDelaunay, Iterators


# Miscellaneous utility funcs

typealias VectorPair{T} NTuple{2, Vector{T}}

"Map a function `f` across rows of the array `A`"
maprows(f, A) = mapslices(f, A, 2)

"Get each adjacent pair from a sequence"
adjacentpairs(itr) = partition(itr, 2, 1)

"2D analogue of the cross product"
crossprod2d(v::Vector, w::Vector) = v[1]*w[2] - v[2]*w[1]

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


# Duplicate-finding

""" Get indices of the unique items in `itr`. For non-unique items, returns the
    highest index. """
unique_indices(itr) = indexin(collect(unique(itr)), collect(itr))


# Mapping to the interval (1+eps(), 2-2*eps())

"Map `x` from the interval [`lb`, `ub`] to (1.0, 2.0)"
norm12(x, lb, ub) = (x - lb) ./ (ub - lb) * (max_coord - min_coord) + min_coord

"Map each element of the vector `x` to the interval (1.0, 2.0)"
norm12(x::Vector) = norm12(x, extrema(x)...)

"Map `x` from the interval (1.0, 2.0) to [`lb`, `ub`]"
unnorm(x, lb, ub) = (x - min_coord) / (max_coord - min_coord) .* (ub - lb) + lb

"Map the log of `x` from [`lb`, `ub`] to the interval (1.0, 2.0)"
lognorm12(x, lb, ub) = norm12(log10(x), log10(lb), log10(ub))

"Map the log of the vector `x` to the interval (1.0, 2.0)"
lognorm12(x::Vector) = norm12(log10(x))

"Map the log of `x` from (1.0, 2.0) back to [`lb`, `ub`]"
unlognorm(x, lb, ub) = 10^unnorm(x, log10(lb), log10(ub))


# Tessellations

typealias Tessellation DelaunayTessellation2D

BoundingBox(t::Tessellation) = BoundingBox(min_coord, max_coord, min_coord, max_coord)

function Base.in(xn, yn, tess::Tessellation)
    if !in(xn, yn, BoundingBox(tess))
        return false
    else
        tri = findtriangle(tess, xn, yn)
        return !isexternal(tri)
    end
end

"Make a Delaunay tessellation on the plane using points in `xs` and `ys`"
function get_tessellation(xs::Vector, ys::Vector; uselog=false)
    xn = uselog ? lognorm12(xs) : norm12(xs) 
    yn = uselog ? lognorm12(ys) : norm12(ys)
    tess = DelaunayTessellation(length(xs))
    for (x, y) in zip(xn, yn)
        push!(tess, Point(x, y))
    end

    tess
end

"Find the triangle in the tessellation `tess` containing points `xn` and `yn`"
function findtriangle(t::Tessellation, xn, yn)
    triangle = locate(t, Point(xn, yn))
end

"Get the vertices of triangle `tri`"
getvertices(tri) = map(f -> f(tri), [geta, getb, getc])

"Get barycentric coordinates on the triangle `tri` from cartesian (`xn`, `yn`)"
function barycoords(tri, xn, yn)
    a, b, c = getvertices(tri)
    x = map(getx, [a, b, c])
    y = map(gety, [a, b, c])

    r3 = [x[3], y[3]]
    r = [xn, yn]
    T = [x[1]-x[3]  x[2]-x[3];
         y[1]-y[3]  y[2]-y[3]]
    λ1, λ2 = T \ (r - r3)
    λ = [λ1, λ2, 1 - λ1 - λ2]
end

"Linear interpolation from vertex values `x` and barycentric coordinates `λ`"
lininterp(x, λ) = λ[1]*x[1] + λ[2]*x[2] + λ[3]*x[3]