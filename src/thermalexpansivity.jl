# thermalexpansivity.jl
# Calculate the thermal expansivity αᵥ of an EOS

using Calculus: gradient


"Gradient of the log of an EOS, ∇(ln(*ρ*(P, T))), at pressure `P` and temperature `T`"
function gradientlog end
function gradientlog(eos::EOS, P::Float64, T::Float64)
    gradient(PT -> log(eos(PT...)), [P, T])
end
function gradientlog(s::StitchedEOS, P::Float64, T::Float64)
    i = findfirst(eos -> (P, T) in eos, s.eoses)

    i == 0 ? NaN : gradientlog(s.eoses[i], P, T)
end

"Thermal expansivity α of an EOS at pressure `P` and temperature `T`"
function thermalexpansivity(eos::EOS, P::Float64, T::Float64)
    -1 * (gradientlog(eos, P, T)⋅[0, 1])
end
