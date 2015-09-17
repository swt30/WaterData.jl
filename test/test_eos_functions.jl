using FactCheck
import WaterData

include("$(WaterData.config.testdata)/function_testvalues.jl")
res = function_testvalues


facts("Functional equations of state") do
    context("TFD") do
        # our sample pressures
        pressures = [0.1e6, 1e6, 10e6] # in bar
        pressures .*= 1e5              # 1 bar = 1e5 Pa

        function test_TFD(A, Z, n, anticipated_densities)
            TFD = WaterData.TFD(Z, A, n)
            @fact map(TFD, pressures) --> roughly(anticipated_densities, rtol=0.01)
        end

        function test_TFD(A, Z, anticipated_densities)
            TFD = WaterData.TFD(Z, A)
            @fact map(TFD, pressures) --> roughly(anticipated_densities, rtol=0.01)
        end

        context("with a single element") do
            context("Fe") do
                test_TFD(55.845, 26, [5900, 8130, 15400])
            end
            context("Bi") do
                test_TFD(208.98, 83, [21800, 26200, 41000])
            end
        end

        context("with multiple elements") do
            context("TiO2") do
                test_TFD([47.867, 15.9994], [22, 8], [1., 2.], [3190, 4920, 10400])
            end
            context("PbS") do
                test_TFD([207.2, 32.065], [82, 16], [12800, 17000, 29600])
            end
        end
    end

    context("Vinet") do
        p = res.vinetpars
        vinet = WaterData.Vinet(p[:ρ₀], p[:K₀], p[:dK₀])
        vinetsamples = res.vinetdata
        Ps = vinetsamples[:, 1]
        ρs = vinetsamples[:, 2]

        for (P, ρ) in zip(Ps, ρs)
            @fact log(vinet(P)) --> roughly(log(ρ), rtol=0.01)
        end
    end

    context("BME3") do
        p = res.bme3pars
        bme3 = WaterData.BME3(p[:ρ₀], p[:K₀], p[:dK₀])
        bme3samples = res.bme3data
        Ps = bme3samples[:, 1]
        ρs = bme3samples[:, 2]

        for (P, ρ) in zip(Ps, ρs)
            @fact log(bme3(P)) --> roughly(log(ρ), rtol=0.01)
        end
    end

    context("MGD thermal expansivity") do
        p = res.mgdpars
        ice_density(molar_volume) = WaterData.h2o_molar_mass / molar_volume
        ρ₀ = ice_density(p[:V₀])

        for (T, datatable) in res.mgddata
            bme3 = WaterData.BME3(ρ₀, p[:K₀], p[:dK₀])
            mgd = WaterData.MGDPressureEOS(bme3,
                p[:T₀], p[:θD₀], p[:γ₀], p[:q], p[:n])
            Ps = datatable[:, 1] * 1e9  # GPa -> Pa
            Vs = datatable[:, 2] / 1e6  # cm3 /mol -> m3 /mol
            ρs = map(ice_density, Vs)    # kg / m3

            for (P, ρ) in zip(Ps, ρs)
                @fact log(mgd(P, T)) --> roughly(log(ρ), rtol=0.01)
            end
        end
    end

    context("Choukroun-Grasset low-temperature ice") do
        cgpars = WaterData
        Ps = res.cgdata[:, 1] * 1e6  # MPa -> Pa
        Vs = res.cgdata[:, 2] / 1e3  # cm3 /g -> m3 /kg
        ρs = 1 ./ Vs              # kg / m3
        Ts = res.cgdata[:, 3]        # K

        cgIp = WaterData.get_choukroungrasset_pars(:ice_I)
        cgIIIp = WaterData.get_choukroungrasset_pars(:ice_III)
        cgVp = WaterData.get_choukroungrasset_pars(:ice_V)
        cgVIp = WaterData.get_choukroungrasset_pars(:ice_VI)

        cgI = WaterData.ChoukrounGrasset(cgIp...)
        cgIII = WaterData.ChoukrounGrasset(cgIIIp...)
        cgV = WaterData.ChoukrounGrasset(cgVp...)
        cgVI = WaterData.ChoukrounGrasset(cgVIp...)

        @fact cgI(Ps[1], Ts[1]) --> roughly(ρs[1], rtol=0.01)
        @fact map(cgIII, Ps[2:3], Ts[2:3]) --> roughly(ρs[2:3], rtol=0.01)
        @fact map(cgV, Ps[4:5], Ts[4:5]) --> roughly(ρs[4:5], rtol=0.01)
        @fact cgVI(Ps[6], Ts[6]) --> roughly(ρs[6], rtol=0.01)
    end

    context("EOS inversion") do
        context("Round-trip upon inversion gives the same values") do
            # the Vinet is a function that works by inversion
            f_fe = WaterData.Vinet(8.30e3, 156.2e9, 6.08)  # BME for iron
            ρ_fe = 10000
            P_fe = WaterData.pressure(f_fe, ρ_fe)
            @fact f_fe(P_fe) --> roughly(ρ_fe)

            # so does the BME
            f_h2o = WaterData.BME(1.46e3, 23.7e9,  4.15)  # BME for ice VII
            ρ_h2o = 5000
            P_h2o = WaterData.pressure(f_h2o, ρ_h2o)
            @fact f_h2o(P_h2o) --> roughly(ρ_h2o)
        end
    end
end
