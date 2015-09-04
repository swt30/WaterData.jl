using FactCheck
import WaterData


facts("Constant verification") do
    @fact WaterData.R_h2o --> roughly(WaterData.R / WaterData.h2o_molar_mass, rtol=0.01)
end
