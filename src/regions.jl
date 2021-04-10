struct BoxTruncationRegion <: TruncationRegion
    a::Vector{Float64}
    b::Vector{Float64}
end

intruncationregion(r::BoxTruncationRegion, x::AbstractArray) = all(r.a .<= x) && all(x .<= r.b) 

abstract type EllipticalTruncationRegion <: TruncationRegion end
abstract type PolytopeTruncationRegion <: TruncationRegion end