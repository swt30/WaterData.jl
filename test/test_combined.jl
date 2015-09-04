using FactCheck
import WaterData


# preliminary type definitions

"type definitions for testing piecewise EOS"
module test_combined_resources
import WaterData

"1D test EOS"
type TestEOS1D <: WaterData.FunctionalEOS
    power::Int
end
Base.call(e::TestEOS1D, P::Real) = P^(e.power)

"2D test EOS"
type TestEOS2D <: WaterData.FunctionalEOS
    power::Int
    add::Int
end
Base.call(e::TestEOS2D, P::Real, T::Real) = (P+T)^(e.power) + (e.add)

end # module resources

import test_combined_resources
res = test_combined_resources


# tests

facts("Combined EOS") do
    context("1D piecewise EOS") do
        # three different versions
        identity = res.TestEOS1D(1)
        squared = res.TestEOS1D(2)
        cubed = res.TestEOS1D(3)

        # make piecewise EOS
        eos = WaterData.PressurePiecewiseEOS(
            [identity, squared, cubed], [4, 7, 9, 11])

        context("Get correct EOS from the piecewise EOS") do
            before4 = WaterData.get_single_eos(eos, 2)
            after11 = WaterData.get_single_eos(eos, 12)
            @fact typeof(before4) --> WaterData.OutOfDomainEOS
            @fact typeof(after11) --> WaterData.OutOfDomainEOS
        end

        context("Evaluate the piecewise EOS to the correct values") do
            @fact eos(2) --> isnan
            @fact eos(14) --> isnan
            @fact eos(5) --> 5
            @fact eos(8) --> 8^2
            @fact eos(10) --> 10^3
        end
    end

    context("2D stitched EOS") do
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

        context("Values outside the domain are NaN") do
            @fact eos(-3, -7) --> isnan
            @fact eos(11, 12) --> isnan
        end

        context("Correct values within the domain") do
            @fact eos(2, 2) --> (2+2)^1 + 0  # first EOS
            @fact eos(2, 8) --> (2+8)^2 + 1  # second EOS
            @fact eos(8, 8) --> (8+8)^3 + 2  # third EOS
        end

        context("Correct priority for overlapping EOS") do
            @fact eos(4, 4) --> (4+4)^1 + 0  # first EOS only
            @fact eos(4, 8) --> (4+8)^2 + 1  # second EOS only
            @fact eos(2, 8) --> (2+8)^2 + 1  # ditto
        end
    end
end
