cd(@__DIR__)
using Pkg; Pkg.activate("..")
using Revise
# Pkg.precompile()

using TruncatedDistributions
using Test
using Distributions
using HCubature

include("exampleDists.jl")

function test_distributions(dist_gen)
       for dg in dist_gen
              d, r, properties = dg()
              dt = BasicBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b)
              @test length(d) == properties["length"]
              @test tp(dt) ≈ properties["tp"]
              @test mean(dt) ≈ properties["mean"]
       end
end

test_distributions(distribution_generators2)
# test_distributions(distribution_generators3)
# test_distributions(distribution_generators4)
# test_distributions(distribution_generators5)