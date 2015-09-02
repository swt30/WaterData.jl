using FactCheck
import WaterData


# FactCheck helper functions (usable in all tests)

"Closure to check if something is between `a` and `b`"
function between(a, b)
    @assert a < b
    x -> greater_than(a)(x) && less_than(b)(x)
end

"Is `x` either `$(WaterData.min_coord)` or `$(WaterData.max_coord)`?"
function either_min_or_max_coord(x)
    isapprox(x, WaterData.min_coord) || isapprox(x, WaterData.max_coord)
end

"Placeholder function so that we can test that something doesn't error"
noerror(x) = true


# Tests on util.jl

facts("FactCheck helper functions") do
    context("`between`") do
        @fact between(3,5)(4) --> true
        @fact between(7,9)(7) --> false
        @fact between(2,3)(3) --> false
        @fact between(-3,5)(-4) --> false
        @fact between(-7,0)(2) --> false
        @fact_throws AssertionError between(2,2)
        @fact_throws AssertionError between(3,-2)
    end

    context("`either_min_or_max_coord` function") do
        @fact WaterData.min_coord --> either_min_or_max_coord
        @fact WaterData.max_coord --> either_min_or_max_coord
        @fact 1.5 --> not(either_min_or_max_coord)
    end

    context("`noerror`") do
        @fact 1 + 1 --> noerror
    end
end

facts("Utility functions") do
    context("Small utility functions") do
        # maprows
        A = [1 2;3 4]
        @fact vec(WaterData.maprows(sum, A)) --> [3, 7]
        @fact vec(WaterData.maprows(mean, A)) --> [1.5, 3.5]

        # adjacent pairs
        B = [1, 7, 2, 8]
        @fact collect(WaterData.adjacentpairs(B)) --> [(1,7), (7,2), (2,8)]

        # 2D cross product
        @fact WaterData.crossprod2d([1, 1], [1, 1]) --> 0
        @fact WaterData.crossprod2d([1,-1], [3, 4]) --> 7
        @fact WaterData.crossprod2d([1, 2], [3, 4]) --> -2

        # intersection testing
        c1 = ([1, 5], [8, 5])
        c2 = ([1, 1], [8, 8])
        c3 = ([3, 1], [3, 8])
        c4 = ([8, 5], [8, 8])
        c5 = ([4, 1], [7, 4])
        yeses = [(c1,c2), (c1,c3), (c1,c4), (c2,c3), (c2,c4)]
        noes  = [(c1,c5), (c2,c5), (c3,c4), (c3,c5), (c4,c5)]
        for lines in yeses
            @fact WaterData.intersects(lines...) --> true
        end
        for lines in noes
            @fact WaterData.intersects(lines...) --> false
        end
    end

    context("Duplicate finding") do
        @fact WaterData.unique_indices([1,2,3,4,5]) --> [1,2,3,4,5]
        @fact WaterData.unique_indices([1,1,1,2]) --> [3,4]
        @fact WaterData.unique_indices([1,1,2,3,4,5,4,5]) --> [2,3,4,7,8]
        @fact length(WaterData.unique_indices([1,1,1,1,1])) --> 1
        @fact length(WaterData.unique_indices([1,2,3,1,2,3,1,2])) --> 3

        dups = [1,2,3,2,1,2,3]
        uniq = WaterData.unique_indices(dups)
        @fact dups[uniq] --> [1,2,3]
    end

    context("Bounding boxes") do
        bb = WaterData.BoundingBox(1,2,0,3)
        @fact (1.5, 1.5) in bb --> true
        @fact (0, 0) in bb --> false
        infbox = WaterData.BoundingBox(-2.5, 2.5, -Inf, Inf)
        @fact (1, 4) in infbox --> true
        @fact (-2, 3e9) in infbox --> true
        @fact (-8, 0) in infbox --> false
    end

    context("Polygon inclusion testing") do
        xs = [1, 5, 5, 4, 3, 2, 1]
        ys = [1, 1, 3, 2, 3, 2, 3]
        poly = WaterData.Polygon(xs, ys)

        # outside bounding box
        @fact (0, 0) in poly --> false
        @fact (6, 2) in poly --> false
        # clearly inside, non-pathological
        @fact (1.5, 1.5) in poly --> true
        @fact (3, 2.5) in poly --> true
        @fact (4.5, 1.5) in poly --> true
        # outside, but within bounding box, non-pathological
        @fact (2, 2.5) in poly --> false
        @fact (4, 2.5) in poly --> false
        # outside, pathogical cases
        @fact (0, 1) in poly --> false
        @fact (0, 2) in poly --> false
        @fact (0, 3) in poly --> false
        @fact (6, 1) in poly --> false
        @fact (6, 2) in poly --> false
        @fact (6, 3) in poly --> false
        # inside, pathological cases
        @fact (3, 2) in poly --> true
        @fact (1.5, 2) in poly --> true

        # one more to convince myself it works
        xs = [1, 2, 3, 4, 5, 5, 4, 3, 2, 1]
        ys = [3, 3, 4, 3, 3, 2, 2, 1, 2, 2]
        poly = WaterData.Polygon(xs, ys)

        # outside bounding box
        @fact (-3, -2) in poly --> false 
        # clearly inside, non-pathological
        @fact (2.5, 2.5) in poly --> true
        # outside, but within bounding box, non-pathological
        @fact (1.5, 3.5) in poly --> false
        # outside, pathogical cases
        @fact (1, 4) in poly --> false
        @fact (0, 3) in poly --> false
        # inside, pathological cases
        @fact (3, 2) in poly --> true
        @fact (3, 3) in poly --> true
    end

    context("Normalisation functions") do
        norm12 = WaterData.norm12
        unnorm = WaterData.unnorm
        lognorm12 = WaterData.lognorm12
        unlognorm = WaterData.unlognorm
        min_coord = WaterData.min_coord
        max_coord = WaterData.max_coord

        context("Normalisation to (1,2)") do
            @fact norm12(4, 4, 6) --> min_coord
            @fact norm12(6, 4, 6) --> max_coord
            @fact norm12(5, 4, 6) --> roughly(1.5)
            @fact norm12(0, 4, 6) --> less_than(min_coord)
            @fact norm12(10, 4, 6) --> greater_than(max_coord)

            r = collect(linspace(20, 30))
            @fact all(norm12(r) .≥ min_coord) --> true
            @fact all(norm12(r) .≤ max_coord) --> true

            @fact lognorm12(10^1.5, 10^1, 10^2) --> roughly(1.5)
            @fact lognorm12(5, 10, 100) --> less_than(min_coord)
            @fact lognorm12(500, 10, 100) --> greater_than(max_coord)
            @fact lognorm12(10, 10, 100) --> min_coord
            @fact lognorm12(100, 10, 100) --> max_coord

            r2 = collect(logspace(3, 4))
            @fact all(lognorm12(r2) .≥ min_coord) --> true
            @fact all(lognorm12(r2) .≤ max_coord) --> true
        end

        context("De-normalising from (1,2)") do
            @fact unnorm(min_coord, 40, 60) --> 40
            @fact unnorm(max_coord, 40, 60) --> 60
            @fact unnorm(norm12(50, 40, 60), 40, 60) --> roughly(50)
            @fact unnorm(0.5, 40, 60) --> less_than(40)
            @fact unnorm(2.5, 40, 60) --> greater_than(60)

            @fact unlognorm(min_coord, 10, 100) --> roughly(10)
            @fact unlognorm(max_coord, 10, 100) --> roughly(100)
            @fact unlognorm(lognorm12(50, 10, 100), 10, 100) --> roughly(50)
            @fact unlognorm(0.5, 10, 100) --> less_than(10)
            @fact unlognorm(500, 10, 100) --> greater_than(500)
        end
    end

    context("Tessellation and interpolation") do
        geta = WaterData.geta
        getb = WaterData.getb
        getc = WaterData.getc
        getx = WaterData.getx
        gety = WaterData.gety

        xs = [0, 1, 0, 1]
        ys = [0, 0, 1, 1]
        zs = [0, 1, 1, 2]
        @fact WaterData.get_tessellation(xs, ys) --> noerror
        tess = WaterData.get_tessellation(xs, ys)
        @fact WaterData.findtriangle(tess, 1.2, 1.2) --> noerror
        @fact WaterData.findtriangle(tess, 1.8, 1.8) --> noerror
        tri1 = WaterData.findtriangle(tess, 1.2, 1.2)
        tri2 = WaterData.findtriangle(tess, 1.8, 1.8)
        p1, p2, p3 = geta(tri1), getb(tri1), getc(tri1)
        q1, q2, q3 = geta(tri2), getb(tri2), getc(tri2)

        for f in [getx, gety], p in [p1, p2, p3, q1, q2, q3]
            @fact f(p) --> either_min_or_max_coord
        end

        λ1 = WaterData.barycoords(tri1, 1.2, 1.2)
        @fact sum(λ1) --> roughly(1)
        λ2 = WaterData.barycoords(tri2, 1.8, 1.8)
        @fact sum(λ2) --> roughly(1)

        x1 = [0, 1, 2]
        x2 = [8, 9, 10]
        @fact WaterData.lininterp(x1, λ1) --> between(0, 2)
        @fact WaterData.lininterp(x2, λ2) --> between(8, 10)
    end

    context("Testing for points in a tessellation") do
        # it's a parallelogram
        tess = WaterData.get_tessellation([0., 1, 1, 2], [0., 0, 1, 1])

        # any point in (1,2)x(1,2) should be inside the bounding box
        @fact in(1.5, 1.5, WaterData.BoundingBox(tess)) --> true
        @fact in(0.5, 1.5, WaterData.BoundingBox(tess)) --> false

        # only points actually in the parallelogram should test true here
        @fact in(1.1, 1.9, tess) --> false  # top left corner
        @fact in(1.5, 1.5, tess) --> true   # the middle
        @fact in(1.9, 1.1, tess) --> false  # bottom right corner
        @fact in(1.1, 1.1, tess) --> true   # bottom left corner
        @fact in(1.9, 1.9, tess) --> true   # bottom right corner
    end
end
