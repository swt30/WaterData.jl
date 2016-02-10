# Work with the equation of state data files

using JLD  # for loading files

# Functions to read in the JLD files for quick EOS access
# Generates JLD files if they're not yet present

"Load a file from the data directory"
function loaddata(filename)
    try
        load(joinpath(config.datadir, filename))
    catch SystemError
        error("File does not exist: try running `WaterData.generate_jlds!()` first")
    end
end

"Get the tabular equations of state from storage"
load_tabular_eoses() = loaddata("eos-tabular.jld")
"Get phase boundary information from storage"
load_phase_boundaries() = loaddata("phase-boundaries.jld")
"Get functional equations of state from storage"
load_functional_eoses() = loaddata("eos-functional.jld")
"Get equations of state defined piecewise from storage"
load_piecewise_eoses() = loaddata("eos-piecewise.jld")
"Get the complete stitched equation of state and thermal expansivity from storage"
load_full_eos() = loaddata("eos-full.jld")

function generate_jlds!()
    info("Generating data/phase-boundaries.jld")
    save_phase_boundaries!()
    info("Generating data/eos-functional.jld")
    save_functional_eoses!()
    info("Generating data/eos-tabular.jld")
    save_tabular_eoses!()
    info("Generating data/eos-piecewise.jld")
    save_piecewise_eoses!()
    info("Generating data/eos-full.jld")
    save_full_eos!()

    nothing
end
