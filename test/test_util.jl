using FactCheck
import WaterData


facts("Utility functions") do
    context("Small utility functions") do
        # maprows
        A = [1 2;3 4]
        @fact vec(WaterData.maprows(sum, A)) --> [3, 7]
        @fact vec(WaterData.maprows(mean, A)) --> [1.5, 3.5]

        # adjacent pairs
        B = [1, 7, 2, 8]
        @fact collect(WaterData.adjacentpairs(B)) --> [(1,7), (7,2), (2,8)]

        # 2D cross product
        @fact WaterData.crossprod2d([1, 1], [1, 1]) --> 0
        @fact WaterData.crossprod2d([1,-1], [3, 4]) --> 7
        @fact WaterData.crossprod2d([1, 2], [3, 4]) --> -2

        # isoutside
        @fact WaterData.isoutside(3, 2, 4) --> false
        @fact WaterData.isoutside(2., 3., 4.) --> true
        @fact WaterData.isoutside(3., 3., 4.) --> false

        # intersection testing
        c1 = ([1, 5], [8, 5])
        c2 = ([1, 1], [8, 8])
        c3 = ([3, 1], [3, 8])
        c4 = ([8, 5], [8, 8])
        c5 = ([4, 1], [7, 4])
        yeses = [(c1,c2), (c1,c3), (c1,c4), (c2,c3), (c2,c4)]
        noes  = [(c1,c5), (c2,c5), (c3,c4), (c3,c5), (c4,c5)]
        for lines in yeses
            @fact WaterData.intersects(lines...) --> true
        end
        for lines in noes
            @fact WaterData.intersects(lines...) --> false
        end
    end

    context("Duplicate finding") do
        @fact WaterData.unique_indices([1,2,3,4,5]) --> [1,2,3,4,5]
        @fact WaterData.unique_indices([1,1,1,2]) --> [3,4]
        @fact WaterData.unique_indices([1,1,2,3,4,5,4,5]) --> [2,3,4,7,8]
        @fact length(WaterData.unique_indices([1,1,1,1,1])) --> 1
        @fact length(WaterData.unique_indices([1,2,3,1,2,3,1,2])) --> 3

        dups = [1,2,3,2,1,2,3]
        uniq = WaterData.unique_indices(dups)
        @fact dups[uniq] --> [1,2,3]
    end
end
