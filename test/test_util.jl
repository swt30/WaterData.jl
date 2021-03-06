using Base.Test
import WaterData


@testset "Utility functions" begin
    @testset "Small utility functions" begin
        # maprows
        A = [1 2;3 4]
        @test vec(WaterData.maprows(sum, A)) == [3, 7]
        @test vec(WaterData.maprows(mean, A)) == [1.5, 3.5]

        # adjacent pairs
        B = [1, 7, 2, 8]
        @test collect(WaterData.adjacentpairs(B)) == [(1,7), (7,2), (2,8)]

        # 2D cross product
        @test WaterData.crossprod2d([1, 1], [1, 1]) == 0
        @test WaterData.crossprod2d([1,-1], [3, 4]) == 7
        @test WaterData.crossprod2d([1, 2], [3, 4]) == -2

        # isoutside
        @test !WaterData.isoutside(3, 2, 4)
        @test WaterData.isoutside(2., 3., 4.)
        @test !WaterData.isoutside(3., 3., 4.)

        # intersection testing
        c1 = ([1, 5], [8, 5])
        c2 = ([1, 1], [8, 8])
        c3 = ([3, 1], [3, 8])
        c4 = ([8, 5], [8, 8])
        c5 = ([4, 1], [7, 4])
        yeses = [(c1,c2), (c1,c3), (c1,c4), (c2,c3), (c2,c4)]
        noes  = [(c1,c5), (c2,c5), (c3,c4), (c3,c5), (c4,c5)]
        for lines in yeses
            @test WaterData.intersects(lines...)
        end
        for lines in noes
            @test !WaterData.intersects(lines...)
        end
    end

    @testset "Duplicate finding" begin
        @test WaterData.unique_indices([1,2,3,4,5]) == [1,2,3,4,5]
        @test WaterData.unique_indices([1,1,1,2]) == [3,4]
        @test WaterData.unique_indices([1,1,2,3,4,5,4,5]) == [2,3,4,7,8]
        @test length(WaterData.unique_indices([1,1,1,1,1])) == 1
        @test length(WaterData.unique_indices([1,2,3,1,2,3,1,2])) == 3

        dups = [1,2,3,2,1,2,3]
        uniq = WaterData.unique_indices(dups)
        @test dups[uniq] == [1,2,3]
    end
end
