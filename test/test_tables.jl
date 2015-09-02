using FactCheck
import WaterData


# Tests on tables.jl

facts("Tabular equations of state") do
    # make some basic EOS grids
    gridtable = WaterData.GridEOS([1.,2,3], [1.,2,3], [1. 1 1; 1 1 1; 1 1 1])
    nongridtable = WaterData.UnstructuredEOS([1.,2,2,3], [1.,1,2,2], [1.,2,2,3])
    
    context("Checking if a point is in a table") do
        @fact in(2, 2, gridtable) --> true
        @fact in(2, 2, WaterData.BoundingBox(gridtable)) --> true
        @fact in(3, 0.5, gridtable) --> false
        @fact in(3, 0.5, WaterData.BoundingBox(gridtable)) --> false
        @fact in(2, 1.5, nongridtable) --> true
        @fact in(2.1, 1.9, nongridtable) --> true
        @fact in(1.9, 1.1, nongridtable) --> true
        @fact in(3, 4, nongridtable) --> false  # outside the bounding box
        @fact in(1.1, 1.9, WaterData.BoundingBox(nongridtable)) --> true
        @fact in(1.1, 1.9, nongridtable) --> false # inside the BB but outside the EOS
    end

    context("Log-normalisation based on the range of an entire table") do
        xn, yn = WaterData.lognorm12(2.5, 1.5, gridtable)
        @fact xn --> between(1, 2)
        @fact yn --> between(1, 2)
        P, T = WaterData.unlognorm(xn, yn, gridtable)
        @fact P --> roughly(2.5)
        @fact T --> roughly(1.5)

        xn, yn = WaterData.lognorm12(2.5, 1.5, nongridtable)
        @fact xn --> between(1, 2)
        @fact yn --> between(1, 2)
        P, T = WaterData.unlognorm(xn, yn, nongridtable)
        @fact P --> roughly(2.5)
        @fact T --> roughly(1.5)
    end

    context("Linear interpolation on unstructured data") do
        # interpolating in the region of interest gives appropriate values
        @fact WaterData.lininterp(nongridtable, 1.9, 1.1) --> between(1, 2)
        @fact WaterData.lininterp(nongridtable, 2.1, 1.9) --> between(2, 3)
        # interpolating outside the bounding box gives NaN
        @fact WaterData.lininterp(nongridtable, 0.1, 0.2) --> isnan
        @fact WaterData.lininterp(nongridtable, 4, 5) --> isnan
        # interpolating outside the tessellation but within the BB gives NaN
        @fact WaterData.lininterp(nongridtable, 1.1, 1.9) --> isnan
        @fact WaterData.lininterp(nongridtable, 2.9, 1.1) --> isnan
    end
end