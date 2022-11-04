using TruncatedDistributions
using Test
using Distributions
using HCubature
using PDMats

include("exampleDists.jl")

N_for_random = 10^6

function test_distributions(dist_gen)
       d, r, properties = dg()

       #Test the BasicBoxTruncatedMvNormal
       d_basic = BasicBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b)
       moment_powers = properties["some_moment_to_check"] #specific array of moment powers to check

       @show length(d_basic)
       @show tp(d_basic)
       @show mean(d_basic)
       @show moment(d_basic, moment_powers)

       @test length(d_basic) == properties["length"]
       @test tp(d_basic) ≈ properties["tp"]
       @test mean(d_basic) ≈ properties["mean"]
       #todo - test cov
       @test moment(d_basic, moment_powers) ≈ properties["some_moment_to_check_value"]

       #Test the RecursiveMomentsBoxTruncatedMvNormal class
       d_recursive = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b)
       @test length(d_recursive) == properties["length"]
       @test tp(d_recursive) ≈ properties["tp"]
       @test mean(d_recursive) ≈ properties["mean"]
       #todo - test cov
       @test moment(d_recursive, moment_powers) ≈ properties["some_moment_to_check_value"]

       #Estimate the moments via Monte Carlo and print
       moment_func(x) = prod(x[i]^moment_powers[i] for i in 1:length(d_recursive))
       estimated = mean(moment_func.([rand(d_recursive) for _ in 1:N_for_random]))
       specified = properties["some_moment_to_check_value"]
       @show estimated, specified
end


# test_distribution.(distribution_generators)

d, r, _ = distribution_generators[5]()
dtrunc = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b; max_moment_levels = 4);
μ_gradient(dtrunc, dtrunc.untruncated.μ, dtrunc.untruncated.μ, dtrunc.untruncated.Σ)

# U_gradient(d, d.untruncated.μ, d.untruncated.μ, d.untruncated.Σ)
