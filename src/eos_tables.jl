# Tabular equations of state

using DataFrames  # for loading data in tabular format with headings
using JLD  # for accessing data in .jld files (HDF5-based)
using VoronoiDelaunay  # for making Delaunay meshes on unstructured data
using Dierckx  # for interpolation

export UnstructuredEOS, GridEOS, LineEOS


# EOS types

"Equation of state table"
abstract TabularEOS <: EOS

"Equation of state table defined by unstructured (P, T, ρ) points"
immutable UnstructuredEOS <: TabularEOS
    P::Vector{Float64}
    T::Vector{Float64}
    ρ::Vector{Float64}
    tess::Tessellation

    function UnstructuredEOS(P, T, ρ)
        u = unique_indices(zip(P, T))
        new(P[u], T[u], ρ[u], get_tessellation(P[u], T[u], uselog=true))
    end
end

"Equation of state table defined on a (P, T) grid"
immutable GridEOS <: TabularEOS
    P::Vector{Float64}
    T::Vector{Float64}
    ρ::Matrix{Float64}
    spline::Spline2D

    GridEOS(P, T, ρ) = new(P, T, ρ, Spline2D(P, T, ρ, kx=1, ky=1))
end

"Equation of state table defined on a pressure grid"
immutable LineEOS <: TabularEOS
    P::Vector{Float64}
    ρ::Vector{Float64}
    spline::Spline1D

    LineEOS(P, ρ) = new(P, ρ, Spline1D(P, ρ, k=1))
end


# Testing for inclusion

"Bounding box for a given TabularEOS"
BoundingBox(e::TabularEOS) = BoundingBox(extrema(e.P)..., extrema(e.T)...)

function Base.in(P, T, eos::UnstructuredEOS)
    if !in(P, T, BoundingBox(eos))
        return false
    else
        xn, yn = lognorm12(P, T, eos)
        return in(xn, yn, eos.tess)
    end
end
Base.in(P, T, eos::GridEOS) = in(P, T, BoundingBox(eos))


# Tagging temperature dependence

istempdependent(::UnstructuredEOS) = true
istempdependent(::GridEOS) = true
istempdependent(::LineEOS) = false


# Save EOS data to file

"Load EOS data from tabular format and save to JLD files"
function save_tabular_eoses!()
    jldopen("$(config.datadir)/eos-tabular.jld", "w") do file
        file["sugimura"] = let
            df = readtable("$(config.rawdata)/Sugimura.eos", allowcomments=true, separator=',')
            P = collect(df[:P]) # in GPa
            T = collect(df[:T]) # in K
            ρ = collect(df[:rho]) # in g/cm^3

            P = P * 1e9 # now in Pa
            ρ = ρ * 1e3 # now in kg/m^3

            UnstructuredEOS(P, T, ρ)
        end

        file["feistelwagner"] = let
            df = readtable("$(config.rawdata)/FeistelWagner.eos", allowcomments=true)
            P = collect(df[:pressure]) # in MPa
            T = collect(df[:temperature]) # in K
            ρ = collect(df[:density]) # in kg/m^3

            P = P * 1e6 # now in Pa

            UnstructuredEOS(P, T, ρ)
        end

        file["french"] = let
            df = readtable("$(config.rawdata)/French.eos", allowcomments=true)
            P = collect(df[:pressure]) # in kbar
            T = collect(df[:temperature]) # in K
            ρ = collect(df[:density]) # in g/cm^3

            P = P * 1e8 # now in Pa
            ρ = ρ * 1e3 # now in kg/m^3

            UnstructuredEOS(P, T, ρ)
        end

        file["iapws"] = let
            df = readtable("$(config.rawdata)/IAPWS.eos", allowcomments=true)
            P = collect(df[:P]) # in MPa
            T = collect(df[:T]) # in K
            ρ = collect(df[:rho]) # in kg/m^3

            P = P * 1e6 # now in Pa

            UnstructuredEOS(P, T, ρ)
        end

        file["seager_dft"] = let
            df = readtable("$(config.rawdata)/seager-h2o-dft.eos")
            P = collect(df[:P]) # in Pa
            ρ = collect(df[:rho]) # in kg/m^3

            LineEOS(P, ρ)
        end
    end

    return nothing
end


# Normalisation

"Log-normalise `P` and `T` to the shifted 2D unit square (1.0, 2.0), (1.0, 2.0)"
function lognorm12(P, T, eos::TabularEOS)
    logP = log10(P)
    logT = log10(T)
    xn = norm12(logP, extrema(log10(eos.P))...)
    yn = norm12(logT, extrema(log10(eos.T))...)

    xn, yn
end

"Un-log-normalise `xn` and `yn` back to appropriate pressure and temperature ranges"
function unlognorm(xn, yn, eos::TabularEOS)
    logP = unnorm(xn, extrema(log10(eos.P))...)
    logT = unnorm(yn, extrema(log10(eos.T))...)
    P = 10^logP
    T = 10^logT

    P, T
end


# Interpolation and evaluation of the EOS

"Get the density from the `eos` at a particular `vertex` in its tessellation"
function ρ_atvertex(vertex, eos)
    xn, yn = getx(vertex), gety(vertex)
    P, T = unlognorm(xn, yn, eos)

    searchrange = 1:length(eos.P)
    i = findfirst(i -> (eos.P[i] ≈ P && eos.T[i] ≈ T), searchrange)
    if i == 0
        return NaN
    else
        return eos.ρ[i]
    end
end

""" Linear interpolation from a given triangle `tri` in the tessellation of the
    equation of state `eos` using normalised coordinates (`xn`,`yn`) """
function lininterp(tri, eos, xn, yn)
    vertices = getvertices(tri)
    ρ = map(v -> ρ_atvertex(v, eos), vertices)
    λ = barycoords(tri, xn, yn)
    return lininterp(λ, ρ)
end

"Linear interpolation of an equation of state `eos` at values `P` and `T`"
function lininterp(eos::UnstructuredEOS, P, T)
    if !in(P, T, BoundingBox(eos))
        return NaN
    end
    xn, yn = lognorm12(P, T, eos)
    tri = findtriangle(eos.tess, xn, yn)
    if isexternal(tri)
        return NaN
    else
        return lininterp(tri, eos, xn, yn)
    end
end

""" Put the `eos` on a grid of a given `resolution`"""
function grid(eos::UnstructuredEOS, resolution=config.grid_resolution)
    P = logspace(extrema(log10(eos.P))..., resolution)
    T = logspace(extrema(log10(eos.T))..., resolution)
    @time ρ = Float64[lininterp(eos, P, T) for P in P, T in T]
    GridEOS(P, T, ρ)
end

""" Take a slice of a `GridEOS` at a given temperature, turning it into
    an isothermal `LineEOS`. """
function slice(eos::GridEOS, T)
    Ps = eos.P
    ρs = map(P -> eos(P, T), Ps)
    LineEOS(Ps, ρs)
end

Base.call(eos::UnstructuredEOS, P, T) = lininterp(eos, P, T)
Base.call(eos::GridEOS, P, T) = evaluate(eos.spline, P, T)
Base.call(eos::LineEOS, P) = evaluate(eos.spline, P)
Base.call(eos::LineEOS, P, T) = evaluate(eos.spline, P)
