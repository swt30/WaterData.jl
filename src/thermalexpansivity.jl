# Calculate the thermal expansivity αᵥ of an EOS

import ForwardDiff, Calculus


"Gradient of the log of an EOS, ∇(ln(*ρ*(P, T))), at pressure `P` and temperature `T`"
function gradientlog(eos::EOS, P, T)
    ForwardDiff.gradient(PT -> log(eos(PT[1], PT[2])), [P, T])
end
function gradientlog(s::StitchedEOS, P, T)
    i = findfirst(eos -> (P, T) in eos, s.eoses)

    i == 0 ? [NaN, NaN] : gradientlog(s.eoses[i], P, T)
end

"Do finite differencing for the gradient of EOSes that don't support autodiff"
function _glog_fallback(eos::EOS, P, T)
    Calculus.gradient(PT -> log(eos(PT...)), [P, T])
end

# Delaunay triangulations need floating point filtering
gradientlog(eos::UnstructuredEOS, P, T) = _glog_fallback(eos, P, T)
# MGD pressure EOS uses quadgk integration, which fails to work
gradientlog(eos::MGDPressureEOS, P, T) = _glog_fallback(eos, P, T)
gradientlog(eos::BoundedEOS{MGDPressureEOS}, P, T) = _glog_fallback(eos, P, T)

"Thermal expansivity α of an EOS at pressure `P` and temperature `T`"
function thermalexpansivity(eos::EOS, P, T)
    -1 * (gradientlog(eos, P, T)⋅[0, 1])
end
