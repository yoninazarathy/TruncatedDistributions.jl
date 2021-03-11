module TruncatedDistributions

using Distributions
using HCubature
using LinearAlgebra

export
    BoxTruncatedMvNormal,
    moment,
    compute_moments,
    alpha,
    raw_moment,
    mean,
    cov,
    rand

include("common.jl")
include("dynamic_univariate_fit.jl")
include("numericalMoments.jl")

end

