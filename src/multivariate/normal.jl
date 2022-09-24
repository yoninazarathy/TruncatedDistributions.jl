"""
A Box Truncated normal distribution with a naive implementation and representation in the state of the mean and covariance. Works well for very low dimensions (e.g. 2,3,4).
"""
const BasicBoxTruncatedMvNormal = TruncatedMvDistribution{MvNormal,BoxTruncationRegion,TruncatedMvDistributionSecondOrderState}

function BasicBoxTruncatedMvNormal( μₑ::Vector{Float64},
                                    Σₑ::PDMat, 
                                    a::Vector{Float64}, 
                                    b::Vector{Float64})
    d = MvNormal(μₑ,Σₑ)
    r = BoxTruncationRegion(a,b)
    TruncatedMvDistribution{MvNormal,BoxTruncationRegion,TruncatedMvDistributionSecondOrderState}(d,r)
end

"""
A Box Truncated normal distribution with a recursive moment computation implementation.
"""
const RecursiveMomentsBoxTruncatedMvNormal = TruncatedMvDistribution{MvNormal,BoxTruncationRegion,BoxTruncatedMvNormalRecursiveMomentsState}

function RecursiveMomentsBoxTruncatedMvNormal(  μₑ::Vector{Float64},
                                                Σₑ::PDMat, 
                                                a::Vector{Float64}, 
                                                b::Vector{Float64};
                                                max_moment_levels::Int = 2)
    d = MvNormal(μₑ,Σₑ)
    r = BoxTruncationRegion(a,b)
    s = BoxTruncatedMvNormalRecursiveMomentsState(d,r,max_moment_levels)
    TruncatedMvDistribution{MvNormal,BoxTruncationRegion,BoxTruncatedMvNormalRecursiveMomentsState}(d,r,s)
end

# const EllipsoidTruncatedMvNormal = TruncatedMvDistribution{MvNormal,EllipticalTruncationRegion}
# const PolytopeTruncatedMvNormal = TruncatedMvDistribution{MvNormal,PolytopeTruncationRegion}

