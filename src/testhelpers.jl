# Helper functions used in test routines


"Test helper functions"
module testhelpers

"Closure to check if something is between `a` and `b`"
function between(a, b)
    @assert a < b
    x -> a < x < b
end

"Placeholder function so that we can test that something doesn't error"
noerror(x) = true

end # module testhelpers
