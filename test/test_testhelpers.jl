using FactCheck
import WaterData
import WaterData.testhelpers: between, noerror


facts("Test helper functions - quis custodiet ipsos custodes?") do

    context("`between`") do
        @fact between(3,5)(4) --> true
        @fact between(7,9)(7) --> false
        @fact between(2,3)(3) --> false
        @fact between(-3,5)(-4) --> false
        @fact between(-7,0)(2) --> false
        @fact_throws AssertionError between(2,2)
        @fact_throws AssertionError between(3,-2)
    end

    context("`noerror`") do
        @fact 1 + 1 --> noerror
    end
end
