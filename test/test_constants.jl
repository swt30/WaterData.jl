using Base.Test
import WaterData

@testset "Constant verification" begin
    @test isapprox(WaterData.R_h2o, WaterData.R / WaterData.h2o_molar_mass, rtol=0.01)
end
