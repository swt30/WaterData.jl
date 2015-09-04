using FactCheck
import WaterData


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
