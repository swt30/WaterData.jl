# The heat capacity of water

using Dierckx


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

Base.call(cp::GridHeatCapacity, P, T) = evaluate(cp.spline, P, T)
Base.call(cp::ConstantHeatCapacity, P) = cp.value
Base.call(cp::ConstantHeatCapacity, P, T) = cp.value
Base.call(cp::TFuncHeatCapacity, T) = cp.func(T)
Base.call(cp::TFuncHeatCapacity, P, T) = cp.func(T)
Base.call(cp::PTFuncHeatCapacity, P, T) = cp.func(P, T)
