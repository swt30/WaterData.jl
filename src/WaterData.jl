# Main module file

"Thermophysical properties of water"
module WaterData

# Include all other parts of the module
include("config.jl")
include("constants.jl")
include("util.jl")
include("regions.jl")
include("phaseboundaries.jl")
include("eos.jl")
include("thermalexpansivity.jl")

include("testhelpers.jl")

"Get the tabular equations of state from storage"
load_tabular_eoses() = load("$(config.datadir)/eos-tabular.jld")
"Get phase boundary information from storage"
load_phase_boundaries() = load("$(config.datadir)/phaseboundaries.jld")
"Get functional equations of state from storage"
load_functional_eoses() = load("$(config.datadir)/eos-functional.jld")
"Get equations of state defined piecewise from storage"
load_piecewise_eoses() = load("$(config.datadir)/eos-piecewise.jld")
"Get the complete stitched equation of state and thermal expansivity from storage"
load_full_eos() = load("$(config.datadir)/eos-full.jld")

end # module WaterData
