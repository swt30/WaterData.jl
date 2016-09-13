# Contains the configuration module, holding global settings for the package


"Configuration options"
module config
# Data directories
"The directory in which .jld data files are stored"
const datadir = normpath(joinpath(dirname(@__FILE__), "..", "data"))
"The directory in which raw data are stored in text format"
const rawdata = joinpath(datadir, "raw")
"The directory in which test data are stored"
const testdata = normpath(joinpath(dirname(@__FILE__), "..", "test", "data"))

# Gridded EOS parameters
"Grid resolution for making gridded equations of state"
const grid_resolution = 256
"Minimum pressure used in constructing the full EOS [Pa]"
const Pmin = 1e7
"Maximum pressure used in constructing the full EOS [Pa]"
const Pmax = 1e14
"Minimum temperature used in constructing the full EOS [K]"
const Tmin = 275
"Maximum temperature used in constructing the full EOS [K]"
const Tmax = 25000

end # module config
