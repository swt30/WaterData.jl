using Base.Test
import WaterData


# preliminary type definitions

"type definitions for testing piecewise EOS"
module test_eos_combined_resources
import WaterData

"1D test EOS"
type TestEOS1D <: WaterData.FunctionalEOS
    power::Int
end
(e::TestEOS1D)(P::Real) = P^(e.power)

"2D test EOS"
type TestEOS2D <: WaterData.FunctionalEOS
    power::Int
    add::Int
end
(e::TestEOS2D)(P::Real, T::Real) = (P+T)^(e.power) + (e.add)

end # module resources


# test

@testset "Combined EOS" begin
    res = test_eos_combined_resources

    @testset "1D piecewise EOS" begin
        # three different versions
        identity = res.TestEOS1D(1)
        squared = res.TestEOS1D(2)
        cubed = res.TestEOS1D(3)

        # make piecewise EOS
        eos = WaterData.PressurePiecewiseEOS(
            [identity, squared, cubed], [4, 7, 9, 11])

        @testset "Get correct EOS from the piecewise EOS" begin
            OutOfDomain = WaterData.OutOfDomainEOS

            @test typeof(WaterData.extracteos(eos, 2)) == OutOfDomain
            @test typeof(WaterData.extracteos(eos, 5)) == res.TestEOS1D
            @test typeof(WaterData.extracteos(eos, 12)) == OutOfDomain
        end

        @testset "Evaluate the piecewise EOS to the correct values" begin
            @test_throws DomainError eos(2)
            @test_throws DomainError eos(14)
            @test eos(5) == 5
            @test eos(8) == 8^2
            @test eos(10) == 10^3
        end
    end

    @testset "2D stitched EOS" begin
        # three diferent versions
        identity = res.TestEOS2D(1, 0)
        sqr1 = res.TestEOS2D(2, 1)
        cube2 = res.TestEOS2D(3, 2)

        # bounding boxes that overlap a bit
        BB1 = WaterData.BoundingBox(0, 5, 0, 5)
        BB2 = WaterData.BoundingBox(0, 7, 3, 10)
        BB3 = WaterData.BoundingBox(3, 10, 3, 10)

        # apply the bounds
        e1 = WaterData.BoundedEOS(identity, BB1)
        e2 = WaterData.BoundedEOS(sqr1, BB2)
        e3 = WaterData.BoundedEOS(cube2, BB3)

        # stitch together
        eos = WaterData.StitchedEOS([e1, e2, e3])

        @testset "Values outside the domain are errors" begin
            @test_throws DomainError eos(-3, -7)
            @test_throws DomainError eos(11, 12)
        end

        @testset "Correct values within the domain" begin
            @test eos(2, 2) == (2+2)^1 + 0  # first EOS
            @test eos(2, 8) == (2+8)^2 + 1  # second EOS
            @test eos(8, 8) == (8+8)^3 + 2  # third EOS
        end

        @testset "Correct priority for overlapping EOS" begin
            @test eos(4, 4) == (4+4)^1 + 0  # first EOS only
            @test eos(4, 8) == (4+8)^2 + 1  # second EOS only
            @test eos(2, 8) == (2+8)^2 + 1  # ditto
        end
    end
end
