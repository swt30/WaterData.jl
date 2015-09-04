using FactCheck


# FactCheck setup

FactCheck.setstyle(:default)
FactCheck.clear_results()


# make sure our helper functions work

include("test_testhelpers.jl")


# Test the package

include("test_combined.jl")
include("test_config.jl")
include("test_constants.jl")
include("test_eos.jl")
include("test_functions.jl")
include("test_phaseboundaries.jl")
include("test_regions.jl")
include("test_tables.jl")
include("test_thermalexpansivity.jl")
include("test_util.jl")


# Pass results to the calling process

FactCheck.exitstatus()
