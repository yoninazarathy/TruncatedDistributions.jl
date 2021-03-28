module TruncatedDistributions

using Distributions
using HCubature
using LinearAlgebra
using PDMats

export
    BoxTruncatedMvNormal,
    moment,
    compute_moments,
    alpha,
    raw_moment,
    mean,
    cov,
    rand

include("commonTypes.jl")
include("multivariate/boxTruncatedMvNormal.jl")
include("multivariate/numericalMoments.jl")
include("univariate/dynamic_univariate_fit.jl")

end #module