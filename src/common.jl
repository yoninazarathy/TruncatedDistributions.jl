abstract type TruncationRegion end

struct BoxTruncationRegion{T<:Real} <: TruncationRegion
    a::AbstractVector{T}
    b::AbstractVector{T}
end

struct TruncatedMvNormal{T<:Real,R<:TruncationRegion}
    μₑ::AbstractVector
    Σₑ::AbstractMatrix
    A::TruncationRegion
    # n::Int64
end

function TruncatedMvNormal( μₑ::AbstractVector{T1},
                            Σₑ::AbstractMatrix{T1},
                            a::AbstractVector{T2},
                            b::AbstractVector{T2}) where {T1<:Real,T2<:Real}
    # n = length(μₑ)
    TruncatedMvNormal{T1,BoxTruncationRegion{T2}}(μₑ,Σₑ,BoxTruncationRegion(a,b))
end

function mean(d::TruncatedMvNormal{T1,BoxTruncationRegion{T2}}) where {T1<:Real,T2<:Real}
    return 8
end

function var(d::TruncatedMvNormal)
    return 8
end
