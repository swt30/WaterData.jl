if VERSION >= v"0.5"
    using Base.Test
else
    using BaseTestNext
end

import WaterData


@testset "No config tests necessary" begin
    nothing
end
