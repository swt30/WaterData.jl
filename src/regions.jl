# regions.jl
# Defining regions on the plane, testing for inclusion, tessellations, and interpolation

using VoronoiDelaunay  # for tessellations


# Region types

"A region on the plane"
abstract Region

"A boxed region [`xmin`, `xmax`] ⊗ [`ymin`, `ymax`]"
immutable BoundingBox <: Region
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end

""" A polygon whose corners are given in the vectors `x` and `y`

    For example, Polygon([0, 1, 1, 0], [0, 0, 1, 1]) defines the unit square.
    """
immutable Polygon <: Region
    x::Vector{Float64}
    y::Vector{Float64}
    boundingbox::BoundingBox

    # we precalculate the bounding box since we'll be using it a lot
    function Polygon(x, y)
        if length(x) == length(y)
            new(x, y, BoundingBox(extrema(x)..., extrema(y)...))
        else
            error("Uneven number of x and y points")
        end
    end
end
"Get the corners of a `Polygon`"
corners(p::Polygon) = zip(p.x, p.y)
BoundingBox(p::Polygon) = p.boundingbox


# Inclusion testing

Base.in(xy, r::Region) = in(xy..., r)
Base.in(xy, eos::EOS) = in(xy..., eos)
Base.in(x, y, bb::BoundingBox) = (bb.xmin <= x <= bb.xmax) && (bb.ymin <= y <= bb.ymax)
function Base.in(x, y, p::Polygon)
    if !((x, y) in p.boundingbox)
        return false
    end

    # O(n) polygon inclusion test - this could definitely be faster
    inpoly = false
    lastcx = p.x[end]
    lastcy = p.y[end]
    for (cx, cy) in corners(p)
        if cy < y && lastcy >= y || lastcy < y && cy >= y
            if (cx + (y - cy)/(lastcy - cy)*(lastcx - cx)) < x
                inpoly = !inpoly
            end
        end
        lastcx, lastcy = cx, cy
    end

    return inpoly
end


# Tesselation utils: mapping to the interval (1+eps(), 2-2*eps())

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
lininterp(x, λ) = λ⋅x
