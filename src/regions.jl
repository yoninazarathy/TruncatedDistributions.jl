"""
Truncation in a box between a and b.
"""
struct BoxTruncationRegion <: TruncationRegion
    a::Vector{Float64}
    b::Vector{Float64}
end

intruncationregion(r::BoxTruncationRegion, x::AbstractArray) = all(r.a .<= x) && all(x .<= r.b) 

"""
An Elliptical Truncation region.
"""
struct EllipticalTruncationRegion <: TruncationRegion
    H::PDMat
    h::Vector{Float64}
    c::Float64
end

intruncationregion(r::EllipticalTruncationRegion, x::AbstractArray) = (x-r.h)'*r.H*(x-r.h) <= r.c


abstract type PolytopeTruncationRegion <: TruncationRegion end