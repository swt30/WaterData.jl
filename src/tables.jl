using DataFrames, JLD, VoronoiDelaunay

"Equation of State table"
abstract EOSTable

# type Hull
#     vertices::Matrix{Float64}
# end
# function Hull(x, y)
#     h = pyhull(hcat(x, y), qhull_options="QbB")
#     vertex_indices = h[:vertices] + 1  # for 1-based indexing
#     vertex_locs = hcat(x[vertex_indices], y[vertex_indices])
#     Hull(vertex_locs)
# end

"Equation of state table defined by unstructured (P, T, ρ) points"
type UnstructuredEOSTable <: EOSTable
    P::Vector{Float64}
    T::Vector{Float64}
    ρ::Vector{Float64}
    tess::DelaunayTessellation2D

    function UnstructuredEOSTable(P, T, ρ)
        u = unique_indices(zip(P, T))
        new(P[u], T[u], ρ[u], get_tesselation(P[u], T[u]))
    end
end

"Equation of state table defined on a (P, T) grid"
type GriddedEOSTable <: EOSTable
    P::Vector{Float64}
    T::Vector{Float64}
    ρ::Matrix{Float64}
end

function save_eoses_to_file()
    sugimura = let
        df = readtable("$(config.datadir)/Sugimura.eos", allowcomments=true, separator=',')
        P = collect(df[:P]) # in GPa
        T = collect(df[:T]) # in K
        ρ = collect(df[:rho]) # in kg/m^3

        P = P * 1e9 # now in Pa

        UnstructuredEOSTable(P, T, ρ)
    end

    feistelwagner = let
        df = readtable("$(config.datadir)/FeistelWagner.eos", allowcomments=true)
        P = collect(df[:pressure]) # in MPa
        T = collect(df[:temperature]) # in K
        ρ = collect(df[:density]) # in kg/m^3

        P = P * 1e6

        UnstructuredEOSTable(P, T, ρ)
    end

    french = let
        df = readtable("$(config.datadir)/French.eos", allowcomments=true)
        P = collect(df[:pressure]) # in kbar
        T = collect(df[:temperature]) # in K
        ρ = collect(df[:density]) # in g/cm^3

        P = P * 1e8 # now in Pa
        ρ = ρ * 1e3 # now in kg/m^3

        UnstructuredEOSTable(P, T, ρ)
    end

    iapws = let
        df = readtable("$(config.datadir)/IAPWS.eos", allowcomments=true)
        P = collect(df[:P]) # in MPa
        T = collect(df[:T]) # in K
        ρ = collect(df[:rho]) # in kg/m^3

        P = P * 1e6 # now in Pa

        UnstructuredEOSTable(P, T, ρ)
    end

    eoses = Dict("sugimura"=>sugimura,
             "feistelwagner"=>feistelwagner,
             "french"=>french,
             "iapws"=>iapws)

    save("$(config.datadir)/eoses.jld", eoses)
end

# eoses = load("$(config.datadir)/eoses.jld")
# const iapws = eoses["iapws"]
# const sugimura = eoses["sugimura"]
# const feistelwagner = eoses["feistelwagner"]
# const french = eoses["french"]

"Get the density from the `eos` at a particular `vertex` in its tesselation"
function ρ_atvertex(vertex, eos)
    xn, yn = getx(vertex), gety(vertex)
    P, T = unlognorm(xn, yn, eos)

    searchrange = 1:length(eos.P)
    i = findfirst(i -> (isapprox(eos.P[i], P) 
                        && isapprox(eos.T[i], T)), searchrange)

    eos.ρ[i]
end

""" Linear interpolation from a triangle `tri` in the tesselation of the
    equation of state `eos` using normalised coordinates (`xn`,`yn`) """
function lininterp(tri, eos, xn, yn)
    vertices = getvertices(tri)
    ρ = map(v -> ρ_atvertex(v, eos), vertices)
    λ = barycoords(tri, xn, yn)
    lininterp(ρ, λ)
end

"Linear interpolation of an equation of state `eos`"
function lininterp(eos, P, T)
    if !(minimum(eos.P) <= P <= maximum(eos.P))
        return NaN
    elseif !(minimum(eos.T) <= T <= maximum(eos.T))
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

""" Put the `eos` on a grid of a given `resolution` (default
    $(config.grid_resolution))"""
function grid(eos::UnstructuredEOSTable, resolution=config.grid_resolution)
    P = logspace(extrema(log10(eos.P))..., resolution)
    T = logspace(extrema(log10(eos.T))..., resolution)
    ρ = Float64[lininterp(eos, P, T) for P in P, T in T]
    GriddedEOSTable(P, T, ρ)
end
