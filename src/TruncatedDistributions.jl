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
    raw_moment_from_indices,
    raw_moment_dict,
    mean,
    cov,
    rand,
    tp,
    TruncatedMvDistributionSecondOrderState,
    BasicBoxTruncatedMvNormal,
    RecursiveMomentsBoxTruncatedMvNormal,
    EllipticalTruncationRegion,
    TruncatedMvDistributionState,
    μ_gradient,
    U_gradient

include("commonTypes.jl")
include("regions.jl")
include("commonOperations.jl")
include("commonCompute.jl")
    include("univariate/distributionsPackageExtensions.jl")

    include("multivariate/boxTruncatedMvNormalRecursiveMomentsState.jl")
    include("multivariate/normal.jl")
    include("multivariate/otherThanNormal.jl")
    include("multivariate/boxTruncatedMvNormalRecursiveMoments.jl")

    include("parameterMatching/dynamicUnivariateFit.jl")
    include("parameterMatching/parameter_gradients.jl")

end #module