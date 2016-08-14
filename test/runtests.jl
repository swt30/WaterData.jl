if VERSION >= v"0.5"
    using Base.Test
else
    using BaseTestNext
end

@testset "WaterData tests" begin
    # make sure our helper functions work
    include("test_testhelpers.jl")

    # Test the package
    include("test_config.jl")
    include("test_constants.jl")
    include("test_eos.jl")
    include("test_eos_combined.jl")
    include("test_eos_functions.jl")
    include("test_eos_tables.jl")
    include("test_heatcapacity.jl")
    include("test_phaseboundaries.jl")
    include("test_regions.jl")
    include("test_thermalexpansivity.jl")
    include("test_util.jl")
end
