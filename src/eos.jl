# Generic equation of state types

export EOS, istempdependent


# EOS type

abstract EOS


# Include EOS subsections

include("eos_functions.jl")
include("eos_tables.jl")
include("eos_combined.jl")


# Testing for inclusion

Base.in(xy, eos::EOS) = in(xy..., eos)


# Marking temperature-dependent EOSes

"Does a given EOS have a temperature dependent component?"
function istempdependent(eos::EOS)
    # by default, this throws an error - later types should override
    # in future versions of Julia, this might be easier to implement as a trait
    throw(NotImplementedError())
end
