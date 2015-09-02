using Calculus, JLD


# Calculating the thermal expansivity of an EOS

function gradientlog(eos::EOS, P::Float64, T::Float64)
    gradient(PT -> log(eos(PT...)), [P, T])
end

function gradientlog(s::StitchedEOS, P::Float64, T::Float64)
    i = findfirst(eos -> (P, T) in eos, s.eoses)
    
    i == 0 ? NaN : gradientlog(s.eoses[i], P, T)
end 

function thermalexpansivity(eos::EOS, P::Float64, T::Float64)
    -1 * (gradientlog(eos, P, T)⋅[0, 1])
end

