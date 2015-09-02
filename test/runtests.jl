using FactCheck


# FactCheck setup

FactCheck.setstyle(:default)
FactCheck.clear_results()


# Test files

include("test_common.jl")
include("test_util.jl")
include("test_tables.jl")
include("test_phaseboundaries.jl")
include("test_functions.jl")
include("test_thermalexpansivity.jl")


# Pass results to the calling process

FactCheck.exitstatus()