using Base.Test
import WaterData

include("$(WaterData.config.testdata)/function_testvalues.jl")
res = function_testvalues


@testset "Functional equations of state" begin
    @testset "TFD" begin
        # our sample pressures
        pressures = [0.1e6, 1e6, 10e6] # in bar
        pressures .*= 1e5              # 1 bar = 1e5 Pa

        function test_TFD(A, Z, n, anticipated_densities)
            TFD = WaterData.TFD(Z, A, n)
            @test isapprox(map(TFD, pressures), anticipated_densities, rtol=0.01)
        end

        function test_TFD(A, Z, anticipated_densities)
            TFD = WaterData.TFD(Z, A)
            @test isapprox(map(TFD, pressures), anticipated_densities, rtol=0.01)
        end

        @testset "with a single element" begin
            @testset "Fe" begin
                test_TFD(55.845, 26, [5900, 8130, 15400])
            end
            @testset "Bi" begin
                test_TFD(208.98, 83, [21800, 26200, 41000])
            end
        end

        @testset "with multiple elements" begin
            @testset "TiO2" begin
                test_TFD([47.867, 15.9994], [22, 8], [1., 2.], [3190, 4920, 10400])
            end
            @testset "PbS" begin
                test_TFD([207.2, 32.065], [82, 16], [12800, 17000, 29600])
            end
        end
    end

    @testset "Vinet" begin
        p = res.vinetpars
        vinet = WaterData.Vinet(p[:ρ₀], p[:K₀], p[:dK₀])
        vinetsamples = res.vinetdata
        Ps = vinetsamples[:, 1]
        ρs = vinetsamples[:, 2]

        for (P, ρ) in zip(Ps, ρs)
            @test isapprox(log(vinet(P)), log(ρ), rtol=0.01)
        end
        @test !WaterData.istempdependent(vinet)
    end

    @testset "BME3" begin
        p = res.bme3pars
        bme3 = WaterData.BME3(p[:ρ₀], p[:K₀], p[:dK₀])
        bme3samples = res.bme3data
        Ps = bme3samples[:, 1]
        ρs = bme3samples[:, 2]

        for (P, ρ) in zip(Ps, ρs)
            @test isapprox(log(bme3(P)), log(ρ), rtol=0.01)
        end
        @test !WaterData.istempdependent(bme3)
    end

    @testset "BME4" begin
        inv_range = [4000., 60000.]
        p = res.bme4pars
        bme4 = WaterData.BME4(p[:ρ₀], p[:K₀], p[:dK₀], p[:d2K₀], inv_range...)
        bme4samples = res.bme4data
        Ps = bme4samples[:, 1]
        ρs = bme4samples[:, 2]

        for (P, ρ) in zip(Ps, ρs)
            @test isapprox(log(bme4(P)), log(ρ), rtol=0.01)
        end
        @test !WaterData.istempdependent(bme4)
    end

    @testset "MGD thermal expansivity" begin
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
                @test isapprox(log(mgd(P, T)), log(ρ), rtol=0.01)
            end
            @test WaterData.istempdependent(mgd)
        end
    end

    @testset "Choukroun-Grasset low-temperature ice" begin
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

        @test isapprox(cgI(Ps[1], Ts[1]), ρs[1], rtol=0.01)
        @test isapprox(map(cgIII, Ps[2:3], Ts[2:3]), ρs[2:3], rtol=0.01)
        @test isapprox(map(cgV, Ps[4:5], Ts[4:5]), ρs[4:5], rtol=0.01)
        @test isapprox(cgVI(Ps[6], Ts[6]), ρs[6], rtol=0.01)
        @test WaterData.istempdependent(cgI)
        @test WaterData.istempdependent(cgIII)
        @test WaterData.istempdependent(cgV)
        @test WaterData.istempdependent(cgVI)
    end

    @testset "EOS inversion" begin
        @testset "Round-trip upon inversion gives the same values" begin
            # the Vinet is a function that works by inversion
            f_fe = WaterData.Vinet(8.30e3, 156.2e9, 6.08)  # BME for iron
            ρ_fe = 10000
            P_fe = WaterData.pressure(f_fe, ρ_fe)
            @test f_fe(P_fe) ≈ ρ_fe
            @test !WaterData.istempdependent(f_fe)

            # so does the BME
            f_h2o = WaterData.BME(1.46e3, 23.7e9,  4.15)  # BME for ice VII
            ρ_h2o = 5000
            P_h2o = WaterData.pressure(f_h2o, ρ_h2o)
            @test f_h2o(P_h2o) ≈ ρ_h2o
            @test !WaterData.istempdependent(f_h2o)
        end
    end
end
