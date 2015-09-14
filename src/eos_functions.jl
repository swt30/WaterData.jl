# Functional equations of state

using DataFrames, JLD, Roots

export ChoukrounGrasset, PolytropicEOS,
    BME, BME3, BME4, Vinet, TFD,
    IAPWS,
    BoundedEOS, MGDPressureEOS


# Shared variables

"Default range of ρ (in kg/m^3) to consider when doing numerical inversion of ρ(P)"
const inversion_ρ_range = [1e-3, 1e8]


# Functional EOS type hierarchy

"An EOS which calculates ρ by evaluating a function ρ = ρ(P[, T])"
abstract FunctionalEOS <: EOS
"An EOS which calculates ρ by numerically inverting a function P = P(ρ[, T])"
abstract InverseFunctionalEOS <: FunctionalEOS


# Specific EOSes

"Choukroun and Grasset's cold ice EOS"
immutable ChoukrounGrasset <: FunctionalEOS
    Pref::Float64
    Tref::Float64
    V₀::Float64
    aT::Vector{Float64}
    P₀::Float64
    aP::Vector{Float64}
end
function ChoukrounGrasset(Pref, Tref, V₀, aT₁, aT₂, aT₃, aT₄, P₀, aP₁, aP₂, aP₃, aP₄)
    ChoukrounGrasset(Pref, Tref, V₀, [aT₁, aT₂, aT₃, aT₄], P₀, [aP₁, aP₂, aP₃, aP₄])
end

function get_choukroungrasset_pars(phase)
    readtable("$(config.rawdata)/ChoukrounGrasset-parameters.dat")[phase]
end

"The polytropic EOS, ρ(P) = P₀ + aP^2"
immutable PolytropicEOS <: FunctionalEOS
    ρ₀::Float64
    a::Float64
    n::Float64
end

"The Birch-Murnaghan EOS"
abstract BME <: InverseFunctionalEOS

"Third-order Birch-Murnaghan EOS"
immutable BME3 <: BME
    ρ₀::Float64
    K₀::Float64
    dK₀::Float64
    ρmin::Float64
    ρmax::Float64
end
BME3(ρ₀, K₀, dK₀) = BME3(ρ₀, K₀, dK₀, inversion_ρ_range...)
BME(ρ₀, K₀, dK₀) = BME3(ρ₀, K₀, dK₀)

"Fourth-order Birch-Murnaghan EOS"
immutable BME4 <: BME
    bme3::BME3
    d2K₀::Float64
    ρmin::Float64
    ρmax::Float64
end
BME4(b::BME3, d2K₀) = BME4(b, d2K₀, inversion_ρ_range...)
BME4(ρ₀, K₀, dK₀, d2K₀) = BME4(BME3(ρ₀, K₀, dK₀), d2K₀)
BME(ρ₀, K₀, dK₀, d2K₀) = BME4(ρ₀, K₀, dK₀, d2K₀)

"The Vinet EOS"
immutable Vinet <: InverseFunctionalEOS
    ρ₀::Float64
    K₀::Float64
    dK₀::Float64
    ρmin::Float64
    ρmax::Float64
end
Vinet(ρ₀, K₀, dK₀) = Vinet(ρ₀, K₀, dK₀, inversion_ρ_range...)

"The Thomas-Fermi-Dirac EOS"
immutable TFD <: FunctionalEOS
    Z::Vector{Int}
    A::Vector{Float64}
    n::Vector{Float64}

    function TFD(Z::Vector{Int}, A::Vector{Float64}, n::Vector{Float64})
        @assert length(Z) == length(A) == length(n)
        new(Z, A, n)
    end
end
# default constructor
function TFD{T<:Real}(Z::Vector{Int}, A::Vector{T}, n=ones(length(Z)))
    TFD(Z, Float64[A], n)
end
# put in arrays if necessary
function TFD(Z::Integer, A)
    TFD(Int[Z], [A])
end

"The IAPWS EOS, functional formulation"
immutable IAPWS <: InverseFunctionalEOS
    c::Vector{Float64}
    d::Vector{Float64}
    t::Vector{Float64}
    n::Vector{Float64}
    α::Vector{Float64}
    β::Vector{Float64}
    γ::Vector{Float64}
    ϵ::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    B::Vector{Float64}
    C::Vector{Float64}
    D::Vector{Float64}
    A::Vector{Float64}
    ρmin::Float64
    ρmax::Float64
end
function IAPWS(data::Matrix{Float64}, ρmin, ρmax)
    @assert size(data) == (56, 14)
    IAPWS([data[:, i] for i=1:14]..., ρmin, ρmax)
end


"Wrapper for indicating that some `FunctionalEOS` is bounded"
immutable BoundedEOS{F<:FunctionalEOS} <: EOS
    eos::F
    extent::Region
end
BoundingBox(b::BoundedEOS) = BoundingBox(b.extent)

"Wrapper for including thermal pressure with a BME or MGD EOS"
immutable MGDPressureEOS <: InverseFunctionalEOS
    eos::Union(BME, Vinet)
    T₀::Float64
    θD₀::Float64
    γ₀::Float64
    q::Int
    n::Int
end


# Helper functions for EOS evaluation

# supporting functions for IAPWS density/temperature
"Get density from dimensionless density"
ρ_from_δ(δ) = δ * ρc
"Get dimensionless density from density"
δ_from_ρ(ρ) = ρ / ρc
"Get temperature from dimensionless temperature"
T_from_τ(τ) = Tc / τ
"Get dimensionless temperature from temperature"
τ_from_T(T) = Tc / T

# supporting functions for Helmholtz energy
# residual parts (not used in this code but kept here for testing purposes)
function ϕr1(I::IAPWS, δ, τ, i)
    let n = I.n[i], d = I.d[i], t = I.t[i]

        n * δ^d * τ^t
    end
end
function ϕr2(I::IAPWS, δ, τ, i)
    let n = I.n[i], d = I.d[i], t = I.t[i], c = I.c[i]

        n * δ^d * τ^t * exp(-δ^c)
    end
end
function ϕr3(I::IAPWS, δ, τ, i)
    let n = I.n[i], d = I.d[i], t = I.t[i], α = I.α[i],
        ϵ = I.ϵ[i], β = I.β[i], γ = I.γ[i]

        n * δ^d * τ^t * exp(-α*(δ-ϵ)^2 - β*(τ-γ)^2)
    end
end
function ϕr4(I::IAPWS, δ, τ, i)
    let n = I.n[i], Δ = Δ(I, δ, τ, i), b = I.b[i], ψ = ψ(I, δ, τ, i)

        n * Δ^b * δ * ψ
    end
end

# helper functions for residual parts
function ψ(I::IAPWS, δ, τ, i)
    let C = I.C[i], D = I.D[i]
        exp(-C*(δ-1)^2 - D*(τ-1)^2)
    end
end
function θ(I::IAPWS, δ, τ, i)
    let A = I.A[i], β = I.β[i]
        (1-τ) + A*((δ-1)^2)^(1/2β)
    end
end
function Δ(I::IAPWS, δ, τ, i)
    let θ = θ(I, δ, τ, i), B = I.B[i], a = I.a[i]

        θ^2 + B*((δ-1)^2)^a
    end
end

# derivatives of residual parts
function dϕr1(I::IAPWS, δ, τ, i)
    let n = I.n[i], d = I.d[i], t = I.t[i]

        n * d * δ^(d-1) * τ^t
    end
end
function dϕr2(I::IAPWS, δ, τ, i)
    let n = I.n[i], d = I.d[i], t = I.t[i], c = I.c[i]

        n * exp(-δ^c) * (δ^(d-1) * τ^t * (d - c*δ^c))
    end
end
function dϕr3(I::IAPWS, δ, τ, i)
    let n = I.n[i], d = I.d[i], t = I.t[i], α = I.α[i],
        ϵ = I.ϵ[i], β = I.β[i], γ = I.γ[i]

        (n * δ^d * τ^t * exp(-α*(δ-ϵ)^2 - β*(τ-γ)^2)
            * (d/δ - 2α*(δ - ϵ)))
    end
end
function dϕr4(I::IAPWS, δ, τ, i)
    let n = I.n[i], Δ = Δ(I, δ, τ, i), b = I.b[i], ψ = ψ(I, δ, τ, i),
        dψdδ = dψdδ(I, δ, τ, i), dΔbdδ = dΔbdδ(I, δ, τ, i)

        n * (Δ^b * (ψ + δ*dψdδ) + dΔbdδ * δ * ψ)
    end
end

# derivatives of helper functions
function dψdδ(I::IAPWS, δ, τ, i)
    let C = I.C[i], ψ = ψ(I, δ, τ, i)

        -2C * (δ - 1) * ψ
    end
end
function dΔbdδ(I::IAPWS, δ, τ, i)
    let b = I.b[i], Δ = Δ(I, δ, τ, i), dΔdδ = dΔdδ(I, δ, τ, i)

        b * Δ^(b - 1) * dΔdδ
    end
end
function dΔdδ(I::IAPWS, δ, τ, i)
    let A = I.A[i], θ = θ(I, δ, τ, i), β = I.β[i], a = I.a[i], B = I.B[i]

        ((δ - 1) * (A * θ * 2/β * ((δ - 1)^2)^(1/2β - 1)
                    + 2B * a * ((δ - 1)^2)^(a - 1)))
    end
end

"IAPWS Helmholtz energy"
function ϕr(I::IAPWS, δ, τ)
    # we take sums over the functions dϕ1, dϕ2...
    part1 = mapreduce(i -> ϕr1(I, δ, τ, i), (+), 0, 1:7)
    part2 = mapreduce(i -> ϕr2(I, δ, τ, i), (+), 0, 8:51)
    part3 = mapreduce(i -> ϕr3(I, δ, τ, i), (+), 0, 52:54)
    part4 = mapreduce(i -> ϕr4(I, δ, τ, i), (+), 0, 55:56)

    # and then add them all together
    (part1 + part2 + part3 + part4)
end

"Derivative of the IAPWS Helmholtz energy in the density direction"
function dϕrδ(I::IAPWS, δ, τ)
    # we take sums over the functions dϕ1, dϕ2...
    part1 = mapreduce(i -> dϕr1(I, δ, τ, i), (+), 0, 1:7)
    part2 = mapreduce(i -> dϕr2(I, δ, τ, i), (+), 0, 8:51)
    part3 = mapreduce(i -> dϕr3(I, δ, τ, i), (+), 0, 52:54)
    part4 = mapreduce(i -> dϕr4(I, δ, τ, i), (+), 0, 55:56)

    # and then add them all together
    (part1 + part2 + part3 + part4)
end

"Pressure component of the Choukroun-Gras￼set ice EOS"
function ξP(cg::ChoukrounGrasset, P)
    let a = cg.aP, P₀ = cg.P₀
        a[1] + a[2]*(P - P₀) + a[3]*(P - P₀)^a[4]
    end
end

"Temperature component of the Choukroun-Grasset ice EOS"
function ξT(cg::ChoukrounGrasset, T)
    let a = cg.aT, Tref = cg.Tref
        a[1] + a[2]*(T - Tref) + a[3]*(T - Tref)^a[4]
    end
end

"Specific volume of the Choukroun-Grasset ice EOS at a given pressure and temperature"
specificvolume(cg::ChoukrounGrasset, P, T) = cg.V₀ * ξP(cg, P) * ξT(cg, T)

"Integrand in the MGD expression"
mgd_integrand(t) = t^3 / (exp(t) - 1)

"Vibrational energy in the MGD expression [J / mol]"
mgd_Evib(T, n, θD) = 9n*R*T*(T/θD)^3 * quadgk(mgd_integrand, 0, θD/T)[1]

"Thermal pressure component of the MGD"
function thermalpressure(mgd::MGDPressureEOS, ρ, T)
    let γ₀ = mgd.γ₀,
        θD₀ = mgd.θD₀,
        ρ₀ = mgd.eos.ρ₀,
        q = mgd.q,
        M = h2o_molar_mass,
        n = mgd.n,
        T₀ = mgd.T₀

        γ = γ₀ * (ρ/ρ₀)^(-q)
        θD = θD₀ * exp((γ₀ - γ)/q)

        ΔP = γ * ρ / M * (mgd_Evib(T, n, θD) - mgd_Evib(T₀, n, θD))
    end
end

const _tfd_g_coeffs = [0          0          0         0         0;
                       0          0          0         0         0;
                       1.512E-2   8.955E-2   1.090E-1  5.089     -5.980;
                       2.181E-3   -4.015E-1  1.698     -9.566    9.873;
                       -3.328E-4  5.167E-1   -2.369    1.349E1   -1.427E1;
                       -1.384E-2  -6.520E-1  3.529     -2.095E1  2.264E1]

"Helper function for TFD partial sum"
function β_(n::Integer, ϵ::Vector{Float64})
    n += 1 # adjust n from 2-5 to 3-6
    let g = _tfd_g_coeffs
        return 1./((g[n, 1]
                  + g[n, 2] * (ϵ.^(1/2))
                  + g[n, 3] * ϵ
                  + g[n, 4] * (ϵ.^(3/2))
                  + g[n, 5] * (ϵ.^2)).^(n-1))
    end
end


# Evaluating the EOSes

"Get the pressure for an EOS: P = pressure(ρ[, T])"
function pressure(b::BME3, ρ)
    let ρ₀ = b.ρ₀, K₀ = b.K₀, dK₀ = b.dK₀, η = ρ/ρ₀
        return  (3/2*K₀*(η^(7/3) - η^(5/3))
                 * (1 + 3/4*(dK₀ - 4)*(η^(2/3) - 1)))
    end
end

function pressure(b::BME4, ρ)
    let ρ₀ = b.bme3.ρ₀, K₀ = b.bme3.K₀, dK₀ = b.bme3.dK₀, d2K₀ = b.d2K₀, η = ρ/ρ₀
        return pressure(b) + (3/2*K₀*(η^(7/3) - η^(5/3))
                              * 3/8*(η^(2/3) - 1)^2
                              * (K₀*d2K₀ + dK₀*(dK₀ - 7) + 143/9))
    end
end

function pressure(v::Vinet, ρ)
    let ρ₀ = v.ρ₀, K₀ = v.K₀, dK₀ = v.dK₀, η = ρ/ρ₀
        return (3K₀*η^(2/3) * (1 - η^(-1/3))
                * exp(3/2*(dK₀ - 1)*(1 - η^(-1/3))))
    end
end

"Get the base pressure (no thermal component) for an EOS with MGD thermal expansion"
basepressure(mgd::MGDPressureEOS, ρ) = pressure(mgd.eos, ρ)
pressure(mgd::MGDPressureEOS, ρ, T) = basepressure(mgd, ρ) + thermalpressure(mgd, ρ, T)

function pressure(I::IAPWS, ρ, T)
    δ = δ_from_ρ(ρ)
    τ = τ_from_T(T)

    P = ρ * R_h2o * T * (1 + δ*dϕrδ(I, δ, τ))
end

pressure(beos::BoundedEOS, ρ) = pressure(beos.eos, ρ)
pressure(beos::BoundedEOS, ρ, T) = pressure(beos.eos, ρ, T)

function Base.call(eos::TFD, P)
    let Z = eos.Z, A = eos.A, n = eos.n
        # P is in Pa but we want it in dyn/cm^2: 1 Pa = 10 dyn/cm^2
        P *= 10

        # pre-calculations
        ζ   = (P / 9.524E13)^(1/5) .* Z.^(-2/3)
        ϵ   = (3 ./(32π^2 .* Z.^2)).^(1/3)
        ϕ   = (3^(1/3))/20 + ϵ./(4 .*(3 .^(1/3)))
        α   = 1./(1.941E-2 - ϵ.^(1/2).*6.277E-2 + ϵ.*1.076)
        x₀0 = (8.884E-3 + (ϵ.^(1/2)).*4.998E-1
                            + ϵ.*5.2604E-1).^(-1)
        β₀ = x₀0.*ϕ - 1
        β₁ = β₀.*α + ((1 + β₀)./ϕ)
        β₂ = β_(2, ϵ)
        β₃ = β_(3, ϵ)
        β₄ = β_(4, ϵ)
        β₅ = β_(5, ϵ)

        βζ = β₀ + β₁.*ζ + β₂.*ζ.^2 + β₃.*ζ.^3 + β₄.*ζ.^4 + β₅.*ζ.^5

        x₀ = 1./(ζ + ϕ) .* (1 + exp(-α.*ζ).*βζ)

        num = sum(n.*A)
        denom = sum(n.*x₀.^3 ./ Z)
        ρ = num/denom * 3.866

        # rho is in g/cm3 but we want it in kg/m3: 1 g/cm3 = 1000 kg/m3
        return 1000ρ
    end
end

Base.call(eos::TFD, P, T) = eos(P)

function Base.call(cg::ChoukrounGrasset, P, T)
    P = P/1e6  # Pa -> MPa
    density = 1./specificvolume(cg, P, T)
end

function Base.call(eos::InverseFunctionalEOS, P)
    fzero(ρ -> pressure(eos, ρ) - P, ρmin(eos), ρmax(eos))
end

function Base.call(eos::InverseFunctionalEOS, P, T)
    fzero(ρ -> pressure(eos, ρ, T) - P, ρmin(eos), ρmax(eos))
end

"Minimum density of an EOS"
ρmin(eos::EOS) = eos.ρmin
ρmin(mgd::MGDPressureEOS) = ρmin(mgd.eos)
"Maximum density of an EOS"
ρmax(eos::EOS) = eos.ρmax
ρmax(mgd::MGDPressureEOS) = ρmax(mgd.eos)

function Base.call(eos::PolytropicEOS, P)
    let a = eos.a, ρ₀ = eos.ρ₀, n = eos.n
        return ρ₀ + a*P^n
    end
end

Base.call(eos::PolytropicEOS, P, T) = eos(P)
Base.call(beos::BoundedEOS, P) = beos.eos(P)
Base.call(beos::BoundedEOS, P, T) = beos.eos(P, T)


# Defining the regions within which each function holds

# Bounded EOSes are defined by a bounding box or polygon region
Base.in(x, y, beos::BoundedEOS) = in(x, y, beos.extent)

# Other EOSes have no limit as to where they may be evaluated
Base.in(x, y, eos::EOS) = true


# Save EOS data to file

"Load EOS data from tabular format and save to JLD files"
function save_functional_eoses!()
    jldopen("$(config.datadir)/eos-functional.jld", "w") do file
        # phase boundaries
        phaseregions = let
            pb = PhaseBoundary

            I = map(phase -> pb(:I, phase), [:L, :III, :II])
            II = map(phase -> pb(:II, phase), [:I, :III, :V, :VI])
            III = map(phase -> pb(:III, phase), [:L, :V, :II, :I])
            V = map(phase -> pb(:V, phase), [:L, :VI, :II, :III])
            VI = map(phase -> pb(:VI, phase), [:L, :VII, :VIII, :II, :V])
            VII = map(phase -> pb(:VII, phase), [:L, :X, :VIII, :VI])
            VIII = map(phase -> pb(:VIII, phase), [:VI, :VII, :X])
            X = map(phase -> pb(:X, phase), [:L, :VIII, :VII])

            # fix up some of those columns so the polygons aren't weird
            function flipcolumns!(phase, cols)
                for c in cols
                    reverse!(phase[c].P)
                    reverse!(phase[c].T)
                end
            end
            flipcolumns!(I, 3)
            flipcolumns!(III, (3, 4))
            flipcolumns!(V, (2, 3))
            flipcolumns!(VI, (2, 3, 4))
            flipcolumns!(VII, 3)
            flipcolumns!(X, (2, 3))

            # join the boundaries together
            concatP(phase) = vcat([pb.P for pb in phase]...)
            concatT(phase) = vcat([pb.T for pb in phase]...)

            regions = Dict(
                 "I" => Polygon(concatP(I), concatT(I)),
                 "II" => Polygon(concatP(II), concatT(II)),
                 "III" => Polygon(concatP(III), concatT(III)),
                 "V" => Polygon(concatP(V), concatT(V)),
                 "VI" => Polygon(concatP(VI), concatT(VI)),
                 "VII" => Polygon(concatP(VII), concatT(VII)),
                 "VIII" => Polygon(concatP(VIII), concatT(VIII)),
                 "X" => Polygon(concatP(X), concatT(X)))
        end

        # Choukroun and Grasset's equations of state
        ckg = let
            cgtable = readtable("$(config.rawdata)/ChoukrounGrasset-parameters.dat")
            regionI = phaseregions["I"]
            regionIII = phaseregions["III"]
            regionV = phaseregions["V"]
            regionVI = phaseregions["VI"]

            Dict(
                "I" => BoundedEOS(ChoukrounGrasset(cgtable[:ice_I]...), regionI),
                "III" => BoundedEOS(ChoukrounGrasset(cgtable[:ice_III]...), regionIII),
                "V" => BoundedEOS(ChoukrounGrasset(cgtable[:ice_V]...), regionV),
                "VI" => BoundedEOS(ChoukrounGrasset(cgtable[:ice_VI]...), regionVI))
        end
        write(file, "choukroungrasset", ckg)

        # Polytropic EOS
        pt = Dict(
            "mgsio3" => PolytropicEOS(4100., 0.00161, 0.541),
            "fe" => PolytropicEOS(8300., 0.00349, 0.528),
            "h2o" => PolytropicEOS(1460., 0.00311, 0.513),
            "graphite" => PolytropicEOS(2250., 0.00350, 0.514),
            "sic" => PolytropicEOS(3220., 0.00172, 0.537))
        write(file, "polytropic", pt)

        # Seager's EOS
        sg = Dict(
            "fe_eps" =>   Vinet(8.30e3, 156.2e9, 6.08),
            "h2o" =>        BME(1.46e3, 23.7e9,  4.15),
            "mgsio3_pv" =>  BME(4.10e3, 247.0e9, 3.97, -0.016e-9),
            "fe_tfd" =>     TFD(26, 55.845),
            "h2o_tfd" =>    TFD([1, 8], [1.00794, 15.9994], [2., 1.]),
            "mgsio3_tfd" => TFD([12, 14, 8], [24.305, 28.0855, 15.9994],
                                [1., 1., 3.]))
        write(file, "seager", sg)

        # Miscellaneous
        misc = let
            # TFD in the ice X region
            iceXbound = phaseregions["X"]
            beyondXbound = BoundingBox(5e10, Inf, -Inf, Inf)
            TFD_iceX = BoundedEOS(sg["h2o_tfd"], iceXbound)
            TFD_beyond = BoundedEOS(sg["h2o_tfd"], beyondXbound)

            # IAPWS extensions
            iapwstable = readdlm("$(config.rawdata)/IAPWS-coeffs.dat", ',', Float64, skipstart=1)
            iapwscoeffs = iapwstable[:, 2:end]
            extent1 = BoundingBox(1e9, 1e12, 273.16, Tc)
            ρrange1 = (1e-6, 1e6)
            extent2 = BoundingBox(0.05e6, 1e9, 1273., 25000.)
            ρrange2 = (1e-6, 1e3)
            extent3 = BoundingBox(1e9, 1e12, Tc, 25000)
            ρrange3 = (1e-6, 1e5)
            iapws_highpressure = BoundedEOS(IAPWS(iapwscoeffs, ρrange1...), extent1)
            iapws_hightemp = BoundedEOS(IAPWS(iapwscoeffs, ρrange2...), extent2)
            iapws_highprestemp = BoundedEOS(IAPWS(iapwscoeffs, ρrange3...), extent3)

            # BME with thermal pressure for ice VII
            ρ₀ = 1464.7   # kg/m3
            K₀ = 23.9e9   # Pa
            T₀ = 300      # K
            dK₀ = 4.2     # dimensionless
            γ₀ = 1.2      # dimensionless
            θD₀ = 1470    # K
            q = -2        # dimensionless
            n = 3         # dimensionless

            bme = BME(ρ₀, K₀, dK₀)
            mgd_bme = MGDPressureEOS(bme, T₀, θD₀, γ₀, q, n)
            bounded_mgd_bme = BoundedEOS(mgd_bme, phaseregions["VII"])

            Dict(
                "tfd_iceX" => TFD_iceX,
                "tfd_beyond" => TFD_beyond,
                "iapws_hightemp" => iapws_hightemp,
                "iapws_highpressure" => iapws_highpressure,
                "iapws_highprestemp" => iapws_highprestemp,
                "mgd_iceVII" => bounded_mgd_bme)
        end
        write(file, "misc", misc)
    end
end
