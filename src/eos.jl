# Equation of state types

export EOS


# EOS type

abstract EOS

# Testing for inclusion
Base.in(xy, eos::EOS) = in(xy..., eos)

include("eos_functions.jl")
include("eos_tables.jl")
include("eos_combined.jl")
