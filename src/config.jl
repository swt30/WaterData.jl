module config

"The directory in which data files (tables etc) are stored"
const datadir = joinpath(dirname(@__FILE__), "..", "data")
const rawdata = joinpath(datadir, "raw")
"Grid resolution for making gridded equations of state"
const grid_resolution = 100

end # module config