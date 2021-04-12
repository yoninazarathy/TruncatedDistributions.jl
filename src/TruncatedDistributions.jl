module TruncatedDistributions

using Distributions
using HCubature
using LinearAlgebra
using PDMats

import Distributions: insupport, pdf, moment
import Base: size, length, show, rand
import Statistics: mean, cov

export 
    compute_tp,
    compute_mean,
    compute_cov,
    compute_moment,
    TruncationRegion,
    TruncatedMvDistribution,
    BoxTruncatedMvNormalRecursiveMoments,#QQQQ
    BoxTruncationRegion,
    intruncationregion,
    insupport,
    BoxTruncatedMvNormal,
    moment,
    moments,
    compute_moments,
    alpha,
    pdf,
    raw_moment,
    mean,
    cov,
    rand,
    tp,
    TruncatedMvDistributionSecondOrderState,
    BasicBoxTruncatedMvNormal,
    RecursiveMomentsBoxTruncatedMvNormal

include("commonTypes.jl")
include("regions.jl")
include("commonOperations.jl")
include("commonCompute.jl")
include("univariate/distributionsPackageExtensions.jl")
include("multivariate/boxTruncatedMvNormalRecursiveMomentsState.jl")
include("multivariate/normal.jl")
include("multivariate/otherThanNormal.jl")
include("parameterMatching/dynamicUnivariateFit.jl")

end #module