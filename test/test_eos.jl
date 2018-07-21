using Base.Test
import WaterData

module test_eos_resources
import WaterData

mutable struct NoTempInfo <: WaterData.EOS; end
notempinfo = NoTempInfo()

end


@testset "EOS base tests" begin
    res = test_eos_resources
    @testset "Tagging whether an EOS is temperature dependent" begin
        @test_throws WaterData.NotImplementedError WaterData.istempdependent(res.notempinfo)
    end
end
