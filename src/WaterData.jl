# Main module file

__precompile__()

"Thermophysical properties of water"
module WaterData

# Include all other parts of the module
include("config.jl")
include("constants.jl")
include("util.jl")
include("regions.jl")
include("phaseboundaries.jl")
include("eos.jl")
include("heatcapacity.jl")
include("thermalexpansivity.jl")
include("io.jl")

include("testhelpers.jl")

end # module WaterData
