using FactCheck
import WaterData

module test_eos_resources
import WaterData

type NoTempInfo <: WaterData.EOS; end
notempinfo = NoTempInfo()

end


facts("EOS base tests") do
    res = test_eos_resources
    context("Tagging whether an EOS is temperature dependent") do
        @fact_throws WaterData.istempdependent(res.notempinfo)
    end
end
