module WaterData

include("common.jl")
include("config.jl")
include("util.jl")
include("constants.jl")
include("phaseboundaries.jl")
include("tables.jl")
include("functions.jl")
include("thermalexpansivity.jl")
include("combined.jl")

load_tablular_eoses() = load("$(config.datadir)/eos-tabular.jld")
load_phaseboundaries() = load("$(config.datadir)/phaseboundaries.jld")
load_functional_eoses() = load("$(config.datadir)/eos-functional.jld")
load_full_eos() = load("$(config.datadir)/eos-full.jld")

end # module
