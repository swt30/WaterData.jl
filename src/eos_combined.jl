# Defines equations of state which are combinations of others.

using JLD  # for accessing data files in .jld files (HDF5-based)
using ProgressMeter  # for showing progress bars during calculations

export PressurePiecewiseEOS, StitchedEOS


# This EOS is used for signalling we're outside the domain
struct OutOfDomainEOS <: EOS end
(::OutOfDomainEOS)(args...) = throw(DomainError())


# Piecewise one-dimensional EOS

"Equation of state which is stored piecewise in some coordinate"
abstract type PiecewiseEOS <: EOS end

""" Equation of state which is stored piecewise in the pressure coordinate

    * `eoses`: Vector of `EOS`es. Any EOS can be chosen, though it's assumed you
      will build a `PressurePiecewiseEOS` from one-dimensional `EOS`es.
    * `edges`: Vector of pressure values that bracket each individual EOS. For
      example, the array [0, 1, 2] would serve to define an EOS with a domain
      [0, 2], piecewise from [0, 1] and (1, 2].

    Calling the `PressurePiecewiseEOS` will evaluate the correct EOS at a given
    pressure. """
struct PressurePiecewiseEOS <: PiecewiseEOS
    eoses::Vector{EOS}
    edges::Vector{Float64}

    function PressurePiecewiseEOS(eoses, edges)
        @assert(length(edges) == length(eoses) + 1,
                "EOS/edge mismatch: need one more edge than number of EOSes")
        new(eoses, edges)
    end
end

"Find the appropriate individual EOS in a 1D piecewise `eos` at point `x`"
function extracteos(eos::PiecewiseEOS, x::Number)
    if isoutside(x, first(eos.edges), last(eos.edges))
        if x ≈ first(eos.edges) || x ≈ last(eos.edges)
            # might be caused by floating point problems
            warn("""Attempted to evaluate a piecewise EOS just outside its domain.
                    This may be a floating point precision error.
                    Consider widening the domain slightly.""")
            end
        return OutOfDomainEOS()
    end

    for (layer_num, layer_end) in enumerate(eos.edges[2:end])
        if x <= layer_end
            # it's within the half-open interval (x, y]
            return eos.eoses[layer_num]
        end
    end
end
extracteos(eos::EOS, args...) = eos  # fallback for when it's not piecewise

# calling a PressurePiecewiseEOS just evaluates the appropriate EOS
(eos::PressurePiecewiseEOS)(P) = extracteos(eos, P)(P)
(eos::PressurePiecewiseEOS)(P, T) = extracteos(eos, P)(P, T)


""" Save piecewise EOSes to `eos-piecewise.jld`:

    * The water EOS used by Sara Seager in her 2007 paper, consisting of an ice
    VII BME at low pressures, an intermediate density functional theory
    calculation, and the TFD at high pressures.
    * A copy of the above which is bounded to the region of (P,T) space defined
    by ice X, continuing to higher pressures.
    * The MgSiO3 EOS from that paper, consisting of a 4th-order BME at low
    pressures and the TFD at high pressures.
    * The Fe EOS from that paper, consisting of a Vinet fit at low pressures and
    the TFD at high pressures."""
function save_piecewise_eoses!()
    funcs = load("$(config.datadir)/eos-functional.jld")
    tables = load("$(config.datadir)/eos-tabular.jld")

    phaseregions = load_phase_boundaries()["regions"]

    "Seager's H2O piecewise EOS"
    h2o = let
        eoses = [funcs["seager"]["h2o"],
                 tables["seager_dft"],
                 funcs["seager"]["h2o_tfd"]]
        transitionP = [0, 44.3e9, 7686e9, Inf]
        raw = PressurePiecewiseEOS(eoses, transitionP)
        P = logspace(5, 19)
        ρ = map(raw, P)
        LineEOS(P, ρ)
    end

    "Seager's MgSiO3 piecewise EOS"
    mgsio3 = let
        eoses = [funcs["seager"]["mgsio3_pv"],
                 funcs["seager"]["mgsio3_tfd"]]
        transitionP = [0, 1.35e13, Inf]
        raw = PressurePiecewiseEOS(eoses, transitionP)
        P = logspace(5, 19)
        ρ = map(raw, P)
        LineEOS(P, ρ)
    end

    "Seager's Fe piecewise EOS"
    fe = let
        eoses = [funcs["seager"]["fe_eps"],
                 funcs["seager"]["fe_tfd"]]
        transitionP = [0, 2.09e14, Inf]
        raw = PressurePiecewiseEOS(eoses, transitionP)
        P = logspace(5, 19)
        ρ = map(raw, P)
        LineEOS(P, ρ)
    end

    iceXbound = phaseregions["X"]
    beyondXbound = BoundingBox(4e10, Inf, -Inf, Inf)
    iceX = BoundedEOS(h2o, iceXbound)
    iceX_beyond = BoundedEOS(h2o, beyondXbound)

    jldopen("$(config.datadir)/eos-piecewise.jld", "w") do file
        addrequire(file, WaterData)
        write(file, "h2o", h2o)
        write(file, "mgsio3", mgsio3)
        write(file, "fe", fe)
        write(file, "iceX", iceX)
        write(file, "iceX_beyond", iceX_beyond)
    end
end


# Stitched (piecewise 2D) EOS

"Equation of state which consists of several other EOS stitched together"
mutable struct StitchedEOS <: EOS
    eoses::Vector{EOS}
end
StitchedEOS(a::EOS, b...) = StitchedEOS([a, b...])

"Bounding box for this EOS"
function BoundingBox(s::StitchedEOS)
    bboxes = [BoundingBox(eos) for eos in s.eoses]
    xmin = minimum([bb.xmin for bb in bboxes])
    xmax = maximum([bb.xmax for bb in bboxes])
    ymin = minimum([bb.ymin for bb in bboxes])
    ymax = maximum([bb.ymax for bb in bboxes])
    BoundingBox(xmin, xmax, ymin, ymax)
end

"Find the appropriate individual EOS in a 2D piecewise `StitchedEOS` at point (`P`,`T`)"
function extracteos(s::StitchedEOS, P, T)
    i = findfirst(e -> (P, T) in e, s.eoses)
    if i == 0
        # no EOS was found for this value of P and T
        OutOfDomainEOS()
    else
        s.eoses[i]
    end
end

# calling a PressurePiecewiseEOS just evaluates the appropriate EOS
(s::StitchedEOS)(P, T) = extracteos(s, P, T)(P, T)


# Full EOS generation

""" Save the full EOS and thermal expansivity to `eos-full.jld`:

    * `raw`: the full `StitchedEOS`, consisting of all H2O tabular and
    functional EOSes
    * `grid`: `GridEOS` version of the above (much faster to evaluate)
    * `gridPlusIdeal`: switches to ideal gas below a certain threshold
    * `thermexp`: the precomputed thermal expansivity for the above """
function save_full_eos!()
    funcs = load_functional_eoses()
    tables = load_tabular_eoses()
    piecewise = load_piecewise_eoses()

    water_idealgas = funcs["misc"]["ideal_gas"]

    # In case of overlapping domains, EOSes listed earlier here have priority
    eos = WaterData.StitchedEOS(tables["iapws"],
                                tables["sugimura"],
                                tables["french"],
                                funcs["choukroungrasset"]["I"],
                                funcs["choukroungrasset"]["III"],
                                funcs["choukroungrasset"]["V"],
                                funcs["choukroungrasset"]["VI"],
                                funcs["misc"]["mgd_iceVII"],
                                funcs["misc"]["iapws_pastfrench"],
                                piecewise["iceX"],
                                piecewise["iceX_beyond"],
                                funcs["misc"]["iapws_highpressure"],
                                funcs["misc"]["iapws_highprestemp"],
                                funcs["misc"]["iapws_hightemp"],
                                funcs["misc"]["fallback"])

    # set the resolution
    Nx = config.grid_resolution
    Ny = config.grid_resolution

    # set the region to be evaluated
    bb = BoundingBox(config.Pmin, config.Pmax, config.Tmin, config.Tmax)
    Ps = logspace(log10(bb.xmin), log10(bb.xmax), Nx)
    Ts = logspace(log10(bb.ymin), log10(bb.ymax), Ny)

    # compute ρ
    ρs = zeros(Nx, Ny)
    meter = Progress(Nx*Ny, "Gridding ρ...")
    for (j, T) in enumerate(Ts), (i, P) in enumerate(Ps)
        ρs[i, j] = eos(P, T)
        next!(meter)
    end

    # compute α
    αs = zeros(Nx, Ny)
    meter = Progress(Nx*Ny, "Calculating α...")
    for (j, T) in enumerate(Ts), (i, P) in enumerate(Ps)
        αs[i, j] = thermalexpansivity(eos, P, T)
        next!(meter)
    end

    # We throw away values less than zero, as we don't physically expect this
    # behaviour (except perhaps for low-pressure ices which have negative
    # expansivity, but these aren't in the region of phase space we're
    # interested in)
    clamp!(αs, 0, Inf)

    # make EOSes and write to file
    # this EOS uses everything except the ideal gas EOS
    grideos = GridEOS(Ps, Ts, ρs)
    # this EOS is the same, but adds the ideal gas below 100bar
    gridPlusIdeal = PressurePiecewiseEOS([water_idealgas, grideos], [0, 100e5, Inf])
    thermexp = GridEOS(Ps, Ts, αs)
    jldopen("$(WaterData.config.datadir)/eos-full.jld", "w") do file
        write(file, "raw", eos)
        write(file, "grid", grideos)
        write(file, "gridPlusIdeal", gridPlusIdeal)
        write(file, "thermexp", thermexp)
    end
end
