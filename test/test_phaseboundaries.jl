using Base.Test
import WaterData

@testset "Phase boundaries" begin
    @testset "have correct long-short code mappings" begin
        iskeyof(dict) = k -> k in keys(dict)
        isvaluein(dict) = v -> v in values(dict)
        phasedict = WaterData.phase_mappings

        @test "liquid" |> iskeyof(phasedict)
        @test "L" |> isvaluein(phasedict)
        @test "L" |> iskeyof(phasedict)
        @test "ice VII" |> iskeyof(phasedict)
        @test "VII" |> iskeyof(phasedict)
        @test "VII" |> isvaluein(phasedict)
        @test phasedict["liquid"] == "L"
        @test phasedict["ice VII"] == "VII"
    end

    @testset "have correct end values on the boundaries" begin
        pb = WaterData.PhaseBoundary

        maxtemp(pb) = maximum(pb.T)
        mintemp(pb) = minimum(pb.T)
        maxpress(pb) = maximum(pb.P)
        minpress(pb) = minimum(pb.P)

        @test_throws BoundsError pb(:L, :II)
        @test_throws BoundsError pb("L", "II")

        @test isapprox(mintemp(pb(:L, :VII)), 355.0, rtol=0.01)
        @test isapprox(maxpress(pb(:L, :VII)), 400000e5, rtol=0.01)
        @test isapprox(mintemp(pb("L", :VII)), 355.0, rtol=0.01)
        @test isapprox(minpress(pb("L", :VII)), 22160e5, rtol=0.01)
        @test isapprox(mintemp(pb(:L, "VII")), 355.0, rtol=0.01)
        @test isapprox(mintemp(pb("L", "VII")), 355.0, rtol=0.01)
        @test isapprox(maxtemp(pb(:VII, :X)), 1500, rtol=0.01)
    end

    @testset "Testing points on and across the boundary" begin
        Ps = [1,2,3,4,5]
        Ts = [5,4,3,2,1]
        pb = WaterData.OtherPhaseBoundary(Ps, Ts)

        line1 = ([1, 1], [5, 5])
        line2 = ([1, 1], [2, 5])
        line3 = ([0.5, 0.5], [0.5, 6])
        line4 = ([4, 4], [6, 6])
        line5 = ([4, 1], [4, 3])
        line6 = ([2, 3], [3, 3])

        @test WaterData.intersects(line1, pb)
        @test WaterData.intersects(line2, pb)
        @test !WaterData.intersects(line3, pb)
        @test !WaterData.intersects(line4, pb)
        @test WaterData.intersects(line5, pb)
        # we treat "just touches the boundary" as "intersects"
        @test WaterData.intersects(line6, pb)
    end
end
