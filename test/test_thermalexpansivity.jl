using Base.Test
import WaterData


"Resources for testing thermal expansivity"
module test_thermalexpansivity_resources

import WaterData

""" A simple EOS to check thermal expansivity determination.

    ρ(P, T) = P^n * T^m so that αᵥ(P, T) = -m/T """
struct SimpleEOS <: WaterData.EOS
    n::Int
    m::Int
end
(e::SimpleEOS)(P, T) = P^(e.n) * T^(e.m)

end # module test_thermalexpansivity_resources

res = test_thermalexpansivity_resources


@testset "Thermal expansivity" begin
    e1 = res.SimpleEOS(1, 1)
    e2 = res.SimpleEOS(2, 3)

    @testset "Get gradient of an EOS" begin
        @test WaterData.gradientlog(e1, 1., 1.) ≈ [1/1, 1/1]
        @test WaterData.gradientlog(e1, 4., 1.) ≈ [1/4, 1/1]
        @test WaterData.gradientlog(e2, 1., 1.) ≈ [2/1, 3/1]
        @test WaterData.gradientlog(e2, 4., 2.) ≈ [2/4, 3/2]
    end

    @testset "Get directional derivative of an EOS" begin
        @test WaterData.partiallnT(e1, 1., 1.) ≈ 1/1
        @test WaterData.partiallnT(e1, 4., 1.) ≈ 1/1
        @test WaterData.partiallnT(e2, 1., 1.) ≈ 3/1
        @test WaterData.partiallnT(e2, 4., 2.) ≈ 3/2
    end
end
