module TruncatedDistributions

using Distributions

export
    TruncatedMvNormal,
    mean


include("common.jl")
include("dynamic_univariate_fit.jl")

end
