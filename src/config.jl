module config

"The directory in which data files (tables etc) are stored"
const datadir = joinpath(dirname(@__FILE__), "..", "data")
"Grid resolution for making gridded equations of state"
const grid_resolution = 100

end # module config