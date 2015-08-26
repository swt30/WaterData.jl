using WaterData
using FactCheck
FactCheck.setstyle(:compact)
FactCheck.clear_results()

include("test_util.jl")
include("test_tables.jl")
include("test_phaseboundaries.jl")

FactCheck.exitstatus()