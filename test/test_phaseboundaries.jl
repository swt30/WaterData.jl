using FactCheck
import WaterData


facts("Phase boundaries") do
    context("have correct long-short code mappings") do
        iskeyof(dict) = k -> k in keys(dict)
        isvaluein(dict) = v -> v in values(dict)
        phasedict = WaterData.phase_mappings

        @fact "liquid" --> iskeyof(phasedict)
        @fact "L" --> isvaluein(phasedict)
        @fact "L" --> iskeyof(phasedict)
        @fact "ice VII" --> iskeyof(phasedict)
        @fact "VII" --> iskeyof(phasedict)
        @fact "VII" --> isvaluein(phasedict)
        @fact phasedict["liquid"] --> "L"
        @fact phasedict["ice VII"] --> "VII"
    end

    context("have correct end values on the boundaries") do
        pb = WaterData.PhaseBoundary

        maxtemp(pb) = maximum(pb.T)
        mintemp(pb) = minimum(pb.T)
        maxpress(pb) = maximum(pb.P)
        minpress(pb) = minimum(pb.P)

        @fact_throws pb(:L, :II)
        @fact_throws pb("L", "II")

        @fact mintemp(pb(:L, :VII)) --> roughly(355.0, rtol=0.01)
        @fact maxpress(pb(:L, :VII)) --> roughly(400000e5, rtol=0.01)
        @fact mintemp(pb("L", :VII)) --> roughly(355.0, rtol=0.01)
        @fact minpress(pb("L", :VII)) --> roughly(22160e5, rtol=0.01)
        @fact mintemp(pb(:L, "VII")) --> roughly(355.0, rtol=0.01)
        @fact mintemp(pb("L", "VII")) --> roughly(355.0, rtol=0.01)
        @fact maxtemp(pb(:VII, :X)) --> roughly(1500, rtol=0.01)
    end

    context("Testing points on and across the boundary") do
        Ps = [1,2,3,4,5]
        Ts = [5,4,3,2,1]
        pb = WaterData.OtherPhaseBoundary(Ps, Ts)

        line1 = ([1, 1], [5, 5])
        line2 = ([1, 1], [2, 5])
        line3 = ([0.5, 0.5], [0.5, 6])
        line4 = ([4, 4], [6, 6])
        line5 = ([4, 1], [4, 3])
        line6 = ([2, 3], [3, 3])

        @fact WaterData.intersects(line1, pb) --> true
        @fact WaterData.intersects(line2, pb) --> true
        @fact WaterData.intersects(line3, pb) --> false
        @fact WaterData.intersects(line4, pb) --> false
        @fact WaterData.intersects(line5, pb) --> true
        # we treat "just touches the boundary" as "intersects"
        @fact WaterData.intersects(line6, pb) --> true
    end
end
