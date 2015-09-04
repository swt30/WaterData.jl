# testhelpers.jl
# Helper functions used in test routines


"Test helper functions and FactCheck extensions"
module testhelpers

using FactCheck
import WaterData

"Closure to check if something is between `a` and `b`"
function between(a, b)
    @assert a < b
    x -> greater_than(a)(x) && less_than(b)(x)
end

"Placeholder function so that we can test that something doesn't error"
noerror(x) = true

end # module testhelpers
