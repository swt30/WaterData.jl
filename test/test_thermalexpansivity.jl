using FactCheck
import WaterData


"Resources for testing thermal expansivity"
module test_thermalexpansivity_resources

import WaterData

""" A simple EOS to check thermal expansivity determination.

    ρ(P, T) = P^n * T^m so that αᵥ(P, T) = -m/T """
type SimpleEOS <: WaterData.EOS
    n::Int
    m::Int
end
Base.call(e::SimpleEOS, P, T) = P^(e.n) * T^(e.m)

end # module test_thermalexpansivity_resources

res = test_thermalexpansivity_resources


facts("Thermal expansivity") do
    e1 = res.SimpleEOS(1, 1)
    e2 = res.SimpleEOS(2, 3)

    context("Get gradient of an EOS") do
        @fact WaterData.gradientlog(e1, 1., 1.) --> roughly([1/1, 1/1])
        @fact WaterData.gradientlog(e1, 4., 1.) --> roughly([1/4, 1/1])
        @fact WaterData.gradientlog(e2, 1., 1.) --> roughly([2/1, 3/1])
        @fact WaterData.gradientlog(e2, 4., 2.) --> roughly([2/4, 3/2])
    end
end
