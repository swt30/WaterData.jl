""" Get indices of the unique items in `itr`. For non-unique items, returns the
    highest index. """
unique_indices(itr) = indexin(collect(unique(itr)), collect(itr))

# Tesselation utilities

"Map `x` from the interval [`lb`, `ub`] to (1.0, 2.0)"
norm12(x, lb, ub) = (x - lb) ./ (ub - lb) * (max_coord - min_coord) + min_coord

"Map each element of the vector `x` to the interval (1.0, 2.0)"
norm12(x::Vector) = norm12(x, extrema(x)...)

"Map `x` from the interval (1.0, 2.0) to [`lb`, `ub`]"
unnorm(x, lb, ub) = (x - min_coord) / (max_coord - min_coord) .* (ub - lb) + lb

"Log-normalise `P` and `T` to the shifted 2D unit square (1.0, 2.0), (1.0, 2.0)"
function lognorm12(P, T, eos)
    logP = log10(P)
    logT = log10(T)
    x = norm12(logP, extrema(log10(eos.P))...)
    y = norm12(logT, extrema(log10(eos.T))...)

    x, y
end

"Un-log-normalise `xn` and `yn` back to appropriate pressure and temperature ranges"
function unlognorm(xn, yn, eos)
    logP = unnorm(xn, extrema(log10(eos.P))...)
    logT = unnorm(yn, extrema(log10(eos.T))...)
    P = 10^logP
    T = 10^logT

    P, T
end

"Make a Delaunay tesselation on the plane using points in `xs` and `ys`"
function get_tesselation(xs::Vector, ys::Vector)
    xn = norm12(xs)
    yn = norm12(ys)
    tess = DelaunayTessellation(length(xs))
    for (x, y) in zip(xn, yn)
        push!(tess, Point(x, y))
    end

    tess
end

"Find the triangle in the tesselation `tess` containing points `xn` and `yn`"
function findtriangle(tesselation, xn, yn)
    triangle = locate(tesselation, Point(xn, yn))
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
lininterp(x, λ) = λ[1]*ρ[1] + λ[2]*ρ[2] + λ[3]*ρ[3]
