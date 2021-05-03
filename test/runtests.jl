# cd(@__DIR__)
# using Pkg; Pkg.activate("..")
# using Revise
# Pkg.precompile()

using TruncatedDistributions
using Test
using Distributions
using HCubature
using PDMats

include("exampleDists.jl")

function test_distributions(dist_gen)
       for dg in dist_gen
              d, r, properties = dg()
              
              d_basic = BasicBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b)
              @test length(d_basic) == properties["length"]
              @test tp(d_basic) ≈ properties["tp"]
              @test mean(d_basic) ≈ properties["mean"]
              
              d_recursive = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b)
              @test length(d_recursive) == properties["length"]
              @test tp(d_recursive) ≈ properties["tp"]
              @test mean(d_recursive) ≈ properties["mean"]
       end
end

test_distributions(distribution_generators2)
# test_distributions(distribution_generators3)
# test_distributions(distribution_generators4)
# test_distributions(distribution_generators5)

# dg = distribution_generators2[1]
# d, r, _ = dg()
# dt = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b)
# mean(dt)


# H = [1 -0.5; 
#     -0.5 1.]

# h = [.0,-0.5]
# c = 5.0


# r = EllipticalTruncationRegion(PDMat(H),h,c)

# x = [0.1,0.3]
# intruncationregion(r,x)

# points = [rand(Uniform(-5,5),2) for _ in 1:10^5]
# in_points = filter((x)->intruncationregion(r,x),points)

# dist = TruncatedMvDistribution{MvNormal,EllipticalTruncationRegion,TruncatedMvDistributionSecondOrderState}(
#        MvNormal([8.5,0],[1 0; 0 5]),
#        EllipticalTruncationRegion(PDMat(H),h,c))

# rand_points = [rand(dist) for _ in 1:10^3]
# scatter(first.(in_points),last.(in_points),legend=false,xlim=(-5,5),ylim=(-5,5))
# scatter!(first.(rand_points),last.(rand_points))