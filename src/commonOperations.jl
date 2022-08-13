"""
QQQQ
"""
function insupport(d::TruncatedMvDistribution{D,R,S}, x::AbstractArray) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    insupport(d.untruncated,x) && intruncationregion(d.region,x)
end

function rand(d::TruncatedMvDistribution{D,R,S}) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    rand_naive(d)
end

function rand_naive(d::TruncatedMvDistribution{D,R,S}) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    while true
        candidate = rand(d.untruncated)
        intruncationregion(d.region,candidate) && return candidate
    end
end

function length(d::TruncatedMvDistribution{D,R,S}) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    length(d.untruncated)
end

function size(d::TruncatedMvDistribution{D,R,S}) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    size(d.untruncated)
end

function pdf(d::TruncatedMvDistribution{D,R,S},x::AbstractArray; worst_tol = 1e-3) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    d.state.tp_err < worst_tol || compute_tp(d)
    if intruncationregion(d.region,x)
        pdf(d.untruncated, x) / d.state.tp
    else
        zeros(length(d))
    end
end

function mean(d::TruncatedMvDistribution{D,R,S}; worst_tol = 1e-3 ) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    d.state.μ_err < worst_tol || compute_mean(d)
    return d.state.μ
end

function cov(d::TruncatedMvDistribution{D,R,S}; worst_tol = 1e-3 ) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    d.state.Σ_err < worst_tol || compute_cov(d)
    return d.state.Σ
end

function moment(d::TruncatedMvDistribution{D,R,S}, k::Vector{Int} ; worst_tol = 1e-3 ) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    compute_moment(d,k)
end

function tp(d::TruncatedMvDistribution{D,R,S}; worst_tol = 1e-3) where {D <: MultivariateDistribution, R <: TruncationRegion, S <: TruncatedMvDistributionState} 
    d.state.tp_err < worst_tol || compute_tp(d)
    d.state.tp
end