abstract type Truncation end
abstract type BoxTruncation <: Truncation end
abstract type EllipticalTruncation <: Truncation end
abstract type TruncatedMultivariateDistribution{T} <: ContinuousMultivariateDistribution where T<:Truncation end

const BoxTruncatedMultivariateDistribution = TruncatedMultivariateDistribution{BoxTruncation}
const EllipsoidTruncatedMultivariateDistribution = TruncatedMultivariateDistribution{EllipticalTruncation}

abstract type AbstractBoxTruncatedMvNormal <: BoxTruncatedMultivariateDistribution end
abstract type AbstractEllipsoidTruncatedMvNormal <: BoxTruncatedMultivariateDistribution end