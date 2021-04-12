"
A truncation region defines a subset of space to which the distribution is truncated. The basic operation supported is `intruncationregion()`, which returns true if a vector is inside the truncation region.
"
abstract type TruncationRegion end

"
A state abstract object representing computed quantities of a truncated multivariate distribution.

Every concrete subtype should expose at least the following two fields.

- `n::Int` The length of the distribution
- `tp::Float64` The probability of falling in the truncation region for the non-truncated case.
- `tp_err::Float64` An estimate of the absolute error of the probability `tp`.

Other subtypes may expose.

- `μ::Vector{Float64}` The mean vector.
- `μ_err::Float64` An estimate of the relative error for the mean vector.
- `Σ::PDMat` The covariance matrix.
- `Σ_err::Float64` An estimate of the relative error for the covariance matrix.

Further subtypes may expose.
- `moment_dict::Dict{Vector{Int},Float64}` A dictionary mapping multivariate moment vectors to estimate quantities.
- `prob_dict` A dictionary mapping probability vectors to estimated quantities. 
"
abstract type TruncatedMvDistributionState end

mutable struct TruncatedMvDistributionSecondOrderState <: TruncatedMvDistributionState
    n::Int 
    tp::Float64
    μ::Vector{Float64}
    Σ::PDMat
    tp_err::Float64
    μ_err::Float64
    Σ_err::Float64
    TruncatedMvDistributionSecondOrderState(d::MultivariateDistribution) = new(     length(d),
                                                                                    NaN,
                                                                                    Vector{Float64}(undef,0),
                                                                                    PDMat(Array{Float64,2}(I,length(d),length(d))),
                                                                                    Inf, Inf, Inf)
end


"A truncated multi-variate distribution composed of a Multivariate Distribution, Truncation Region and a State object implementing computable state."
struct TruncatedMvDistribution{D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    untruncated::D
    region::R
    state::S
end

function TruncatedMvDistribution{D,R,S}(d::D,r::R) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState}
    TruncatedMvDistribution(d,r,S(d))
end

# function tp(d::TruncatedMvDistribution{D,R,TruncatedMvDistributionSecondOrderState}) 
#     where {D <: MultivariateDistribution, R <: TruncationRegion}
#     (d.state.tp, d.state.tp_err)
# end

# function mean(d::TruncatedMvDistribution{D,R,S}) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState}
#     (d.state.μ, d.state.μ_err)
# end

# function cov(d::TruncatedMvDistribution{D,R,S}) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState}
#     (d.s)
# end

    # function moment(d::TruncatedMvDistribution{D,R,S}, k::Vector{Int}) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState}

# function Base.show(io::IO, d::BoxTruncatedMvNormalRecursiveMomentsState) 
#     println(io, "Box Truncated MvNormal")
#     println(io, "n = $(d.n)")
#     println(io, "μₑ = $(d.μₑ)" )
#     println(io, "Σₑ = $(d.Σₑ)" )
#     println(io, "α = $(alpha(d))")
#     println("Limits:")
    
#     for i in 1:d.n
#         println(io, "$i:\t ",(d.a[i],d.b[i]))
#     end
    
#     if d.momentsComputed
#         println("Moments:")
#         for k in keys(d.momentDict)
#             println(io, k, "\t", moment(d,k))
#         end
#         println(io,"mean:", mean(d))
#         println(io,"cov:", cov(d))
#     else
#         println("Moments not computed")
#     end
# end