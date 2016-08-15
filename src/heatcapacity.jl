# The heat capacity of water

using Dierckx  # for interpolation

export HeatCapacity, GridHeatCapacity


# Types for different heat capacity behaviours

"An isobaric heat capacity (cₚ)"
abstract HeatCapacity
"A heat capacity that's not constant"
abstract VaryingHeatCapacity <: HeatCapacity
"A heat capacity that depends varies according to some functional form"
abstract FunctionalHeatCapacity <: VaryingHeatCapacity

"A constant heat capacity (no T/P variation)"
immutable ConstantHeatCapacity{T<:Real} <: HeatCapacity
    value::T
end

"A heat capacity which is a function of temperature"
immutable TFuncHeatCapacity{F} <: FunctionalHeatCapacity
    func::F
end

"A heat capacity which is a function of pressure and temperature"
immutable PTFuncHeatCapacity{F} <: FunctionalHeatCapacity
    func::F
end

"A heat capacity interpolated from a log-linear grid"
immutable GridHeatCapacity <: HeatCapacity
    P::Vector{Float64}
    T::Vector{Float64}
    spline::Spline2D

    function GridHeatCapacity(P, T, c_p)
        new(P, T, Spline2D(P, T, c_p, kx=1, ky=1))
    end
end


# Making heat capacities from files

"Generate a heat capacity from a file"
function GridHeatCapacity(filename, scale=1000)
    data = readdlm(filename, Float64)
    T = vec(data[2:end, 1])
    P = vec(data[1, 2:end])
    c_p = Matrix{Float64}(data[2:end, 2:end])' * scale  # note the transpose
    # default scaling factor is 1000, corresponding to the values in the file
    # being provided in g/cm^3/K (cgs) and converted to kg/m^3/K (SI base)

    GridHeatCapacity(P, T, c_p)
end


# Evaluating heat capacities

(cp::GridHeatCapacity)(P, T) = evaluate(cp.spline, P, T)
(cp::ConstantHeatCapacity)(P) = cp.value
(cp::ConstantHeatCapacity)(P, T) = cp.value
(cp::TFuncHeatCapacity)(T) = cp.func(T)
(cp::TFuncHeatCapacity)(P, T) = cp.func(T)
(cp::PTFuncHeatCapacity)(P, T) = cp.func(P, T)


# Make and save heat capacity

"Save the water heat capacity to a JLD file"
function save_heat_capacity!()
    jldopen("$(config.datadir)/heatcapacity.jld", "w") do file
        # make the heat capacity
        c_p = GridHeatCapacity("$(config.rawdata)/heatcap-h2o.dat")
        write(file, "heatcap_h2o", c_p)
    end
end
