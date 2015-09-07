# Defines equations of state which are combinations of others.

using JLD, ProgressMeter

export PressurePiecewiseEOS, StitchedEOS


# Piecewise 1D EOS

"Equation of state which is stored piecewise in some coordinate"
abstract PiecewiseEOS <: EOS

""" Equation of state which is stored piecewise in the pressure coordinate

    * `eoses`: Vector of `EOS`es. Any EOS can be chosen, though it's assumed you
    will build a `PressurePiecewiseEOS` from 1D EOSes.
    * `edges`: Vector of pressure values that bracket each individual EOS. For
    example, the list [0, 1, 2] would serve to define an EOS with a domain
    [0, 2], piecewise from [0, 1] and (1, 2].

    Calling the `PressurePiecewiseEOS` will evaluate the correct EOS for a given
    pressure. """
immutable PressurePiecewiseEOS <: PiecewiseEOS
    eoses::Vector{EOS}     # any EOS can go in the middle
    edges::Vector{Float64}

    function PressurePiecewiseEOS(eoses, edges)
        @assert length(edges) == length(eoses) + 1
        new(eoses, edges)
    end
end

""" Equation of state representing an area outside the evaluation domain.
    Returns NaN when called. """
immutable OutOfDomainEOS <: EOS
end

"Find the appropriate individual EOS in a 1D piecewise `eos` at point `x`"
function get_single_eos(eos::PiecewiseEOS, x::Real)
    # if outside the domain, take the edge point
    if x < first(eos.edges) || x > last(eos.edges)
        return OutOfDomainEOS()
    end

    for (layer_num, layer_end) in enumerate(eos.edges[2:end])
        if x <= layer_end
            # it's within the half-open interval (x, y]
            return eos.eoses[layer_num]
        end
    end
end

# calling a PressurePiecewiseEOS just evaluates the appropriate EOS
Base.call(eos::PressurePiecewiseEOS, P::Real) = get_single_eos(eos, P)(P)
# OutOfDomainEOS gives NaN when called
Base.call(o::OutOfDomainEOS, P::Real) = NaN
Base.call(o::OutOfDomainEOS, P::Real, T::Real) = NaN

""" Save piecewise EOSes to `eos-functional.jld`:

    * The water EOS used by Sara Seager in her 2007 paper, consisting of an ice
    VII BME at low pressures, an intermediate density functional theory
    calculation, and the TFD at high pressures """
function save_piecewise_eos!()
    funcs = load("$(config.datadir)/eos-functional.jld")
    tables = load("$(config.datadir)/eos-tabular.jld")

    eoses = [funcs["seager"]["h2o"],
    tables["seager_dft"],
    funcs["seager"]["h2o_tfd"]]
    transition_pressures = [0, 44.3e9, 7686e9, Inf]

    raw = PressurePiecewiseEOS(eoses, transition_pressures)
    P = logspace(5, 14)
    ρ = map(raw, P)
    grid = LineEOS(P, ρ)

    seager = Dict("seager-h2o-raw"=>raw, "seager-h2o-grid"=>grid)
    save("$(config.datadir)/eos-piecewise.jld", seager)
end


# Stitched (piecewise 2D) EOS

"Equation of state which consists of several other EOS stitched together"
type StitchedEOS <: EOS  # TODO: should this be a subtype of PiecewiseEOS instead?
    eoses::Vector{EOS}
end
StitchedEOS(a::EOS, b...) = StitchedEOS([a, b...])

"Bounding box for this EOS"
function BoundingBox(s::StitchedEOS)
    bboxes = [BoundingBox(e) for e in s.eoses]
    xmin = minimum([bb.xmin for bb in bboxes])
    xmax = maximum([bb.xmax for bb in bboxes])
    ymin = minimum([bb.ymin for bb in bboxes])
    ymax = maximum([bb.ymax for bb in bboxes])
    BoundingBox(xmin, xmax, ymin, ymax)
end

"Find the appropriate individual EOS in a 2D piecewise `eos` at point (`P`,`T`)"
function get_single_eos(s::StitchedEOS, P, T)
    i = findfirst(e -> (P, T) in e, s.eoses)
    if i == 0
        OutOfDomainEOS()
    else
        s.eoses[i]
    end
end

# calling a PressurePiecewiseEOS just evaluates the appropriate EOS
Base.call(s::StitchedEOS, P, T) = get_single_eos(s, P, T)(P, T)


# Full EOS

""" Save the full EOS and thermal expansivity to `eos-full.jld`:

    * `raw`: the full `StitchedEOS`, consisting of all H2O tabular and
    functional EOSes
    * `grid`: `GridEOS` version of the above (much faster to evaluate)
    * `thermexp`: the precomputed thermal expansivity for the above """
function save_full_eos!()
    funcs = load("$(config.datadir)/eos-functional.jld")
    tables = load("$(config.datadir)/eos-tabular.jld")

    # EOSes listed earlier here have priority
    eos = WaterData.StitchedEOS(
    tables["iapws"],
    tables["sugimura"],
    tables["french"],
    funcs["choukroungrasset"]["I"],
    funcs["choukroungrasset"]["III"],
    funcs["choukroungrasset"]["V"],
    funcs["choukroungrasset"]["VI"],
    funcs["misc"]["mgd_iceVII"],
    funcs["misc"]["tfd_iceX"],
    funcs["misc"]["tfd_beyond"],
    funcs["misc"]["iapws_highpressure"],
    funcs["misc"]["iapws_highprestemp"],
    funcs["misc"]["iapws_hightemp"])

    # set the resolution
    Nx = config.grid_resolution
    Ny = config.grid_resolution

    # set the region to be evaluated
    bb = BoundingBox(config.Pmin, config.Pmax, config.Tmin, config.Tmax)
    Ps = logspace(log10(bb.xmin), log10(bb.xmax), Nx)
    Ts = logspace(log10(bb.ymin), log10(bb.ymax), Ny)

    # compute ρ and α
    ρs = zeros(Nx, Ny)
    αs = zeros(Nx, Ny)
    meter = Progress(Nx*Ny, "Calculating ρ...")
    for (i, P) in enumerate(Ps), (j, T) in enumerate(Ts)
        ρs[i, j] = eos(P, T)
        next!(meter)
    end
    meter = Progress(Nx*Ny, "Calculating α...")
    for (i, P) in enumerate(Ps), (j, T) in enumerate(Ts)
        αs[i, j] = thermalexpansivity(eos, P, T)
        next!(meter)
    end

    # make EOSes and write to file
    grideos = GridEOS(Ps, Ts, ρs)
    thermexp = GridEOS(Ps, Ts, αs)
    save("$(WaterData.config.datadir)/eos-full.jld",
    "raw", eos, "grid", grideos, "thermexp", thermexp)
end
