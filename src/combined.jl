using JLD, ProgressMeter

# Combined EOS

function save_combined_eos!()
    funcs = load("$(config.datadir)/eos-functional.jld")
    tables = load("$(config.datadir)/eos-tabular.jld")

    eos = WaterData.StitchedEOS(
          funcs["choukroungrasset"]["I"],
          funcs["choukroungrasset"]["III"],
          funcs["choukroungrasset"]["V"],
          funcs["choukroungrasset"]["VI"],
          tables["iapws"],
          funcs["misc"]["mgd_iceVII"],
          tables["french"],
          funcs["misc"]["tfd_iceX"],
          funcs["misc"]["tfd_beyond"],
          funcs["misc"]["iapws_highpressure"],
          funcs["misc"]["iapws_highprestemp"],
          funcs["misc"]["iapws_hightemp"])

    Nx = 200
    Ny = 200

    bb = BoundingBox(10^(4.75), 1e14, 275, 25000)
    Ps = logspace(log10(bb.xmin), log10(bb.xmax), Nx)
    Ts = logspace(log10(bb.ymin), log10(bb.ymax), Ny)
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

    grideos = GridEOS(Ps, Ts, ρs)
    thermexp = GridEOS(Ps, Ts, αs)

    save("$(WaterData.config.datadir)/eos-full.jld", 
         "raw", eos, "grid", grideos, "thermexp", thermexp)
end


