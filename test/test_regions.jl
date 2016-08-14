if VERSION >= v"0.5"
    using Base.Test
else
    using BaseTestNext
end

import WaterData
import WaterData.testhelpers: between, noerror

@testset "Regions and tessellations" begin
    @testset "Bounding boxes" begin
        bb = WaterData.BoundingBox(1,2,0,3)
        @test (1.5, 1.5) ∈ bb
        @test (0, 0) ∉ bb
        infbox = WaterData.BoundingBox(-2.5, 2.5, -Inf, Inf)
        @test (1, 4) ∈ infbox
        @test (-2, 3e9) ∈ infbox
        @test (-8, 0) ∉ infbox
    end

    @testset "Polygon inclusion testing" begin
        xs = [1, 5, 5, 4, 3, 2, 1]
        ys = [1, 1, 3, 2, 3, 2, 3]
        poly = WaterData.Polygon(xs, ys)

        # outside bounding box
        @test (0, 0) ∉ poly
        @test (6, 2) ∉ poly
        # clearly inside, non-pathological
        @test (1.5, 1.5) ∈ poly
        @test (3, 2.5) ∈ poly
        @test (4.5, 1.5) ∈ poly
        # outside, but within bounding box, non-pathological
        @test (2, 2.5) ∉ poly
        @test (4, 2.5) ∉ poly
        # outside, pathogical cases
        @test (0, 1) ∉ poly
        @test (0, 2) ∉ poly
        @test (0, 3) ∉ poly
        @test (6, 1) ∉ poly
        @test (6, 2) ∉ poly
        @test (6, 3) ∉ poly
        # inside, pathological cases
        @test (3, 2) ∈ poly
        @test (1.5, 2) ∈ poly

        # one more to convince myself it works
        xs = [1, 2, 3, 4, 5, 5, 4, 3, 2, 1]
        ys = [3, 3, 4, 3, 3, 2, 2, 1, 2, 2]
        poly = WaterData.Polygon(xs, ys)

        # outside bounding box
        @test (-3, -2) ∉ poly
        # clearly inside, non-pathological
        @test (2.5, 2.5) ∈ poly
        # outside, but within bounding box, non-pathological
        @test (1.5, 3.5) ∉ poly
        # outside, pathogical cases
        @test (1, 4) ∉ poly
        @test (0, 3) ∉ poly
        # inside, pathological cases
        @test (3, 2) ∈ poly
        @test (3, 3) ∈ poly
    end

    @testset "Normalisation functions" begin
        norm12 = WaterData.norm12
        unnorm = WaterData.unnorm
        lognorm12 = WaterData.lognorm12
        unlognorm = WaterData.unlognorm
        min_coord = WaterData.min_coord
        max_coord = WaterData.max_coord

        @testset "Normalisation to (1,2)" begin
            @test norm12(4, 4, 6) == min_coord
            @test norm12(6, 4, 6) == max_coord
            @test norm12(5, 4, 6) ≈ 1.5
            @test norm12(0, 4, 6) < min_coord
            @test norm12(10, 4, 6) > max_coord

            r = collect(linspace(20, 30))
            @test all(norm12(r) .≥ min_coord)
            @test all(norm12(r) .≤ max_coord)

            @test lognorm12(10^1.5, 10^1, 10^2) ≈ 1.5
            @test lognorm12(5, 10, 100) < min_coord
            @test lognorm12(500, 10, 100) > max_coord
            @test lognorm12(10, 10, 100) == min_coord
            @test lognorm12(100, 10, 100) == max_coord

            r2 = collect(logspace(3, 4))
            @test all(lognorm12(r2) .≥ min_coord)
            @test all(lognorm12(r2) .≤ max_coord)
        end

        @testset "De-normalising from (1,2)" begin
            @test unnorm(min_coord, 40, 60) == 40
            @test unnorm(max_coord, 40, 60) == 60
            @test unnorm(norm12(50, 40, 60), 40, 60) ≈ 50
            @test unnorm(0.5, 40, 60) < 40
            @test unnorm(2.5, 40, 60) > 60

            @test unlognorm(min_coord, 10, 100) ≈ 10
            @test unlognorm(max_coord, 10, 100) ≈ 100
            @test unlognorm(lognorm12(50, 10, 100), 10, 100) ≈ 50
            @test unlognorm(0.5, 10, 100) < 10
            @test unlognorm(500, 10, 100) > 500
        end
    end

    @testset "Tessellation and interpolation" begin
        geta = WaterData.geta
        getb = WaterData.getb
        getc = WaterData.getc
        getx = WaterData.getx
        gety = WaterData.gety

        xs = [0, 1, 0, 1]
        ys = [0, 0, 1, 1]
        zs = [0, 1, 1, 2]
        @test WaterData.get_tessellation(xs, ys) |> noerror
        tess = WaterData.get_tessellation(xs, ys)
        @test WaterData.findtriangle(tess, 1.2, 1.2) |> noerror
        @test WaterData.findtriangle(tess, 1.8, 1.8) |> noerror
        tri1 = WaterData.findtriangle(tess, 1.2, 1.2)
        tri2 = WaterData.findtriangle(tess, 1.8, 1.8)
        p1, p2, p3 = geta(tri1), getb(tri1), getc(tri1)
        q1, q2, q3 = geta(tri2), getb(tri2), getc(tri2)

        for f in [getx, gety], p in [p1, p2, p3, q1, q2, q3]
            @test f(p) ∈ (WaterData.min_coord, WaterData.max_coord)
        end

        λ1 = WaterData.barycoords(tri1, 1.2, 1.2)
        @test sum(λ1) ≈ 1
        λ2 = WaterData.barycoords(tri2, 1.8, 1.8)
        @test sum(λ2) ≈ 1

        x1 = [0, 1, 2]
        x2 = [8, 9, 10]
        @test WaterData.lininterp(λ1, x1) |> between(0, 2)
        @test WaterData.lininterp(λ2, x2) |> between(8, 10)
    end

    @testset "Testing for points in a tessellation" begin
        # it's a parallelogram
        tess = WaterData.get_tessellation([0., 1, 1, 2], [0., 0, 1, 1])

        # any point in (1,2)x(1,2) should be inside the bounding box
        @test in(1.5, 1.5, WaterData.BoundingBox(tess))
        @test in(0.5, 1.5, WaterData.BoundingBox(tess)) == false

        # only points actually in the parallelogram should test true here
        @test in(1.1, 1.9, tess) == false  # top left corner
        @test in(1.5, 1.5, tess)    # the middle
        @test in(1.9, 1.1, tess) == false  # bottom right corner
        @test in(1.1, 1.1, tess)    # bottom left corner
        @test in(1.9, 1.9, tess)    # bottom right corner
    end
end
