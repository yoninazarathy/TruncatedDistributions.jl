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