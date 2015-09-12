# Calculate the thermal expansivity αᵥ of an EOS

using Calculus: derivative, gradient


"Directional derivative d(ln ρ)/dT"
function partiallnT(eos::EOS, P, T)
    # ForwardDiff.gradient(PT -> log(eos(PT[1], PT[2])), [P, T])[2]
    derivative(δT -> log(eos(P, T + δT)), 0)::Float64
end
function partiallnT(s::StitchedEOS, P, T)
    i = findfirst(eos -> (P, T) in eos, s.eoses)
    i == 0 ? NaN : partiallogy(s.eoses[i], P, T)
end

"Gradient ∇(ln ρ)"
function gradientlog(eos::EOS, P, T)
    gradient(PT -> log(eos(PT[1], PT[2])), [P, T])::Vector{Float64}
end
function gradientlog(eos::StitchedEOS, P, T)
    i = findfirst(eos -> (P, T) in eos, s.eoses)
    i == 0 ? [NaN, NaN] : gradient(s.eoses[i], P, T)
end

"Thermal expansivity α of an EOS at pressure `P` and temperature `T`"
function thermalexpansivity(eos::EOS, P, T)
    -1 * partiallnT(eos, P, T)
end
