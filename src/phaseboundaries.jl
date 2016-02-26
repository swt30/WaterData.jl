# Calculate and use phase boundaries for water

using Dierckx  # for interpolations

export PhaseBoundary


# Phase boundary names and parameters

"A dictionary mapping phase names to short codes, like 'L' for liquid"
const phase_mappings = let
    keys = ["liquid", "ice I", "ice II", "ice III",
            "ice V", "ice VI", "ice VII", "ice VIII", "ice X"]
    shortkeys = split("L I II III V VI VII VIII X")

    phasemap = Dict(zip(keys, shortkeys))
    shortmap = Dict(zip(shortkeys, shortkeys))

    merge!(phasemap, shortmap)
end

"Read in phase boundary parameter table from Dunaeva et al"
function read_phase_boundary_table()
    readdlm("$(config.rawdata)/Dunaeva-phase-boundaries.dat")
end

"""Describes the T/P extent and shape parameters of a phase boundary"""
immutable PhaseBoundaryPars
    phase1::ASCIIString
    phase2::ASCIIString
    Tmin::Float64
    Tmax::Float64
    Pmin::Float64
    Pmax::Float64
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
end
"Read a phase boundary from a file"
function PhaseBoundaryPars(phase1::AbstractString, phase2::AbstractString)
    table = read_phase_boundary_table()

    match11 = (table[:, 1] .== phase1)
    match22 = (table[:, 2] .== phase2)
    match21 = (table[:, 2] .== phase1)
    match12 = (table[:, 1] .== phase2)
    matches_normal = match11 & match22
    matches_reversed = match12 & match21
    matches = matches_normal | matches_reversed

    row = table[matches, :]
    row[5:6] *= 1e5 # convert from bar to Pa

    PhaseBoundaryPars(row...)
end


# Phase boundary definitions

"Holds details about the boundary between two phases"
abstract PhaseBoundary

"A phase boundary following the formulation of Dunaeva et al"
immutable DunaevaPhaseBoundary <: PhaseBoundary
    P::Vector{Float64}
    T::Vector{Float64}
    pars::PhaseBoundaryPars
end
PhaseBoundary(P::AbstractVector, T, pars) = DunaevaPhaseBoundary(collect(P), T, pars)

"A phase boundary that doesn't follow the Dunaeva formulation"
immutable OtherPhaseBoundary <: PhaseBoundary
    P::Vector{Float64}
    T::Vector{Float64}
    spline::Spline1D

    OtherPhaseBoundary(P, T) = new(P, T, Spline1D(P, T, k=1, bc="error"))
end

"Minimum pressure of the phase boundary (padded to avoid problems at P=0)"
Pmin(pbp::PhaseBoundaryPars) = pbp.Pmin == 0 ? 1e-9 : pbp.Pmin
Pmin(dpb::DunaevaPhaseBoundary) = Pmin(dpb.pars)
Pmin(pb::OtherPhaseBoundary) = minimum(pb.P)

"Maximum pressure of the phase boundary"
Pmax(pbp::PhaseBoundaryPars) = pbp.Pmax == 0 ? 1e-9 : pbp.Pmin
Pmax(dpb::DunaevaPhaseBoundary) = Pmax(dpb.pars)
Pmax(pb::OtherPhaseBoundary) = maximum(pb.P)

"Pressure along a phase boundary"
pressure(pb::PhaseBoundary) = pb.P

getpoint(pb::PhaseBoundary, i) = [pb.P[i], pb.T[i]]
Base.length(pb::PhaseBoundary) = length(pb.P)
Base.start(pb::PhaseBoundary) = 1
Base.next(pb::PhaseBoundary, i) = getpoint(pb, i), i+1
Base.done(pb::PhaseBoundary, i) = i > length(pb)

"Calculate a temperature at a given pressure `P` using Dunaeva parameters `pars`"
function dunaevatemp(pars::PhaseBoundaryPars, P)
    P = P / 100000
    T = pars.a + pars.b*P + pars.c*log(P) + pars.d/P + pars.e*sqrt(P)
end

"Temperature along a phase boundary"
temperature(pb::PhaseBoundary) = pb.T
temperature(pb::OtherPhaseBoundary, P) = pb.spline(P)
temperature(pb::DunaevaPhaseBoundary, P) = dunaevatemp(pb.pars, P)

# Construct a phase boundary between phases `phase1` and `phase2`
function PhaseBoundary(phase1::AbstractString, phase2::AbstractString)
    pars = PhaseBoundaryPars(phase1, phase2)
    P = collect(linspace(pars.Pmin, pars.Pmax))
    T = map(T -> dunaevatemp(pars, T), P)

    PhaseBoundary(P, T, pars)
end
function PhaseBoundary(phase1::Symbol, phase2::Symbol)
    p1 = phase_mappings[string(phase1)]
    p2 = phase_mappings[string(phase2)]
    PhaseBoundary(p1, p2)
end
function PhaseBoundary(phase1::Symbol, phase2::AbstractString)
    p1 = phase_mappings[string(phase1)]
    PhaseBoundary(p1, phase2)
end
function PhaseBoundary(phase1::AbstractString, phase2::Symbol)
    p2 = phase_mappings[string(phase2)]
    PhaseBoundary(phase1, p2)
end

"Does a given line segment ([x₁, y₁], [x₂, y₂]) intersect a given phase boundary?"
function intersects(AB::VectorPair, pb::PhaseBoundary)
    linesegments = adjacentpairs(pb)

    any(CD -> intersects(AB, CD), linesegments)
end


"Save phase boundaries to file"
function save_phase_boundaries!()
    pb = PhaseBoundary

    I = map(phase -> pb(:I, phase), [:L, :III, :II])
    II = map(phase -> pb(:II, phase), [:I, :III, :V, :VI])
    III = map(phase -> pb(:III, phase), [:L, :V, :II, :I])
    V = map(phase -> pb(:V, phase), [:L, :VI, :II, :III])
    VI = map(phase -> pb(:VI, phase), [:L, :VII, :VIII, :II, :V])
    VII = map(phase -> pb(:VII, phase), [:L, :X, :VIII, :VI])
    VIII = map(phase -> pb(:VIII, phase), [:VI, :VII, :X])
    X = map(phase -> pb(:X, phase), [:L, :VIII, :VII])

    # since the regions are just lists of P,T pairs for each side,
    # we need to ensure that they have a consistent orientation;
    # we do this by flipping the sides that have the wrong orientation
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

    # now that the regions have consistent orientations, we can
    # join the boundaries together
    concatP(phase) = vcat([pb.P for pb in phase]...)
    concatT(phase) = vcat([pb.T for pb in phase]...)

    regions = Dict("I" => Polygon(concatP(I), concatT(I)),
                   "II" => Polygon(concatP(II), concatT(II)),
                   "III" => Polygon(concatP(III), concatT(III)),
                   "V" => Polygon(concatP(V), concatT(V)),
                   "VI" => Polygon(concatP(VI), concatT(VI)),
                   "VII" => Polygon(concatP(VII), concatT(VII)),
                   "VIII" => Polygon(concatP(VIII), concatT(VIII)),
                   "X" => Polygon(concatP(X), concatT(X)))

    dunaeva_table = read_phase_boundary_table()
    dunaeva_boundaries = maprows(dunaeva_table) do row
        p1, p2 = row[1:2]
        PhaseBoundary(p1, p2)
    end

    iapws_boundary = let
        table = readdlm("$(config.rawdata)/iapws-phase-boundary.dat")
        P = table[:, 1] * 1e6 # MPa -> Pa
        T = table[:, 2]
        OtherPhaseBoundary(P, T)
    end

    boundaries = Dict("iapws"=>iapws_boundary,
                      "dunaeva"=>dunaeva_boundaries,
                      "regions"=>regions)

    jldopen("$(config.datadir)/phase-boundaries.jld", "w") do file
        write(file, "boundaries", boundaries)
        write(file, "regions", regions)
    end

    return nothing
end
