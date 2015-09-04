"Thermophysical properties of water"
module WaterData

include("config.jl")
include("constants.jl")
include("util.jl")
include("eos.jl")
include("regions.jl")
include("phaseboundaries.jl")
include("tables.jl")
include("functions.jl")
include("combined.jl")
include("thermalexpansivity.jl")

include("testhelpers.jl")

load_tabular_eoses() = load("$(config.datadir)/eos-tabular.jld")
load_phase_boundaries() = load("$(config.datadir)/phaseboundaries.jld")
load_functional_eoses() = load("$(config.datadir)/eos-functional.jld")
load_full_eos() = load("$(config.datadir)/eos-full.jld")
load_piecewise_eoses() = load("$(config.datadir)/eos-piecewise.jld")

end # module
