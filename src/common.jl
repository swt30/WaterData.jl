# parent equation of state type

abstract EOS
Base.in(xy, eos::EOS) = in(xy..., eos)


# combining EOSes

type StitchedEOS <: EOS
    eoses::Vector{EOS}
end
StitchedEOS(a::EOS, b...) = StitchedEOS([a, b...])

function Base.call(s::StitchedEOS, P, T)
    i = findfirst(e -> (P, T) in e, s.eoses)
    if i == 0
        return NaN
    else
        return s.eoses[i](P, T)
    end
end


# Defining regions on the plane and testing for inclusion

abstract Region

immutable BoundingBox <: Region
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end

function BoundingBox(s::StitchedEOS)
    bboxes = [BoundingBox(e) for e in s.eoses]
    xmin = minimum([bb.xmin for bb in bboxes])
    xmax = maximum([bb.xmax for bb in bboxes])
    ymin = minimum([bb.ymin for bb in bboxes])
    ymax = maximum([bb.ymax for bb in bboxes])
    BoundingBox(xmin, xmax, ymin, ymax)
end

immutable Polygon <: Region
    x::Vector{Float64}
    y::Vector{Float64}
    boundingbox::BoundingBox

    function Polygon(x, y)
        if length(x) == length(y)
            new(x, y, BoundingBox(extrema(x)..., extrema(y)...))
        else 
            error("Uneven number of x and y points")
        end
    end
end
corners(p::Polygon) = zip(p.x, p.y)
BoundingBox(p::Polygon) = p.boundingbox

Base.in(x, y, bb::BoundingBox) = (bb.xmin <= x <= bb.xmax) && (bb.ymin <= y <= bb.ymax)
Base.in(xy, bb::Region) = in(xy..., bb)

function Base.in(x, y, p::Polygon)
    if !((x, y) in p.boundingbox)
        return false
    end

    # O(n) polygon inclusion test
    inpoly = false
    lastcx = p.x[end]
    lastcy = p.y[end]
    for (cx, cy) in corners(p)
        if cy < y && lastcy >= y || lastcy < y && cy >= y
            if (cx + (y - cy)/(lastcy - cy)*(lastcx - cx)) < x
                inpoly = !inpoly
            end
        end
        lastcx, lastcy = cx, cy
    end

    return inpoly
end
