if VERSION >= v"0.5"
    using Base.Test
else
    using BaseTestNext
end

import WaterData
import WaterData.testhelpers: between, noerror

# quis custodiet ipsos custodes?
@testset "Test helper functions" begin
    @testset "`between`" begin
        @test between(3,5)(4)
        @test !between(7,9)(7)
        @test !between(2,3)(3)
        @test !between(-3,5)(-4)
        @test !between(-7,0)(2)
        @test_throws AssertionError between(2,2)
        @test_throws AssertionError between(3,-2)
    end

    @testset "`noerror`" begin
        @test noerror(1 + 1)
    end
end
