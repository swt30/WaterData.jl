module WaterData

include("common.jl")
include("config.jl")
include("util.jl")
include("constants.jl")
include("phaseboundaries.jl")
include("tables.jl")
include("functions.jl")

eoses = load("$(config.datadir)/eoses.jld")
const iapws = eoses["iapws"]
const sugimura = eoses["sugimura"]
const feistelwagner = eoses["feistelwagner"]
const french = eoses["french"]

phaseboundaries = load("$(config.datadir)/phase-boundaries.jld")
const evaporation = phaseboundaries["iapws"]
const freezing = phaseboundaries["dunaeva"]

end # module
