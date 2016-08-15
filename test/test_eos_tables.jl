using Base.Test
import WaterData
import WaterData.testhelpers: between


# Tests on tables.jl

@testset "Tabular equations of state" begin
    # make some basic EOS grids
    gridtable = WaterData.GridEOS([1.,2,3], [1.,2,3], [1. 1 1; 1 1 1; 1 1 1])
    nongridtable = WaterData.UnstructuredEOS([1.,2,2,3], [1.,1,2,2], [1.,2,2,3])
    linetable = WaterData.LineEOS([1., 2, 3], [5., 7, 8])

    @testset "Checking if a point is in a table" begin
        @test in(2, 2, gridtable)
        @test in(2, 2, WaterData.BoundingBox(gridtable))
        @test !in(3, 0.5, gridtable)
        @test !in(3, 0.5, WaterData.BoundingBox(gridtable))
        @test in(2, 1.5, nongridtable)
        @test in(2.1, 1.9, nongridtable)
        @test in(1.9, 1.1, nongridtable)
        @test !in(3, 4, nongridtable)  # outside the bounding box
        @test in(1.1, 1.9, WaterData.BoundingBox(nongridtable))
        @test !in(1.1, 1.9, nongridtable) # inside the BB but outside the EOS
    end

    @testset "Log-normalisation based on the range of an entire table" begin
        xn, yn = WaterData.lognorm12(2.5, 1.5, gridtable)
        @test xn |> between(1, 2)
        @test yn |> between(1, 2)
        P, T = WaterData.unlognorm(xn, yn, gridtable)
        @test P ≈ 2.5
        @test T ≈ 1.5

        xn, yn = WaterData.lognorm12(2.5, 1.5, nongridtable)
        @test xn |> between(1, 2)
        @test yn |> between(1, 2)
        P, T = WaterData.unlognorm(xn, yn, nongridtable)
        @test P ≈ 2.5
        @test T ≈ 1.5
    end

    @testset "Linear interpolation on unstructured data" begin
        # interpolating in the region of interest gives appropriate values
        @test WaterData.lininterp(nongridtable, 1.9, 1.1) |> between(1, 2)
        @test WaterData.lininterp(nongridtable, 2.1, 1.9) |> between(2, 3)
        # interpolating outside the bounding box gives NaN
        @test WaterData.lininterp(nongridtable, 0.1, 0.2) |> isnan
        @test WaterData.lininterp(nongridtable, 4, 5) |> isnan
        # interpolating outside the tessellation but within the BB gives NaN
        @test WaterData.lininterp(nongridtable, 1.1, 1.9) |> isnan
        @test WaterData.lininterp(nongridtable, 2.9, 1.1) |> isnan
    end

    @testset "Linear interpolation of 1-D data" begin
        @test linetable(1) == 5
        @test linetable(1.5) == 6
        @test linetable(2) == 7
        @test linetable(2.5) == 7.5
        @test linetable(3) == 8
        @test linetable(-9) == 5  # clips at start
        @test linetable(28) == 8  # clips at end
    end

    @testset "Checking tagging of temperature dependence" begin
        @test WaterData.istempdependent(gridtable)
        @test !WaterData.istempdependent(linetable)
        @test WaterData.istempdependent(nongridtable)
    end

end
