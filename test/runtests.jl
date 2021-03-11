cd(@__DIR__)
using Pkg; Pkg.activate("..")
using Revise
# Pkg.precompile()

using TruncatedDistributions
using Test

# @testset "TruncatedDistributions.jl" begin
#     # Write your tests here.
#     @test 1 == 1
# end

# μₑ = [3.5,2,3.5,3.5,5.3]
# Σₑ = [ 7. 1   0 0 1 ; 
#        1  3.3 2 0 0 ;
#        0 2 3.4 0 0 ;
#        0 0 0 4 0; 
#        1 0 0 0 2]
# a = [-20.,-20,-20,-20,-20]
# b = [20.,20.,20,20,20]

μₑ = [2.5, 3.5 ,2.3]
Σₑ = [2.0 -0.5 0;
       -0.5 5.0 0;
        0 0 2.3]
a = [-10.,-10.,-10.]
b = [10.,10.,10.]


d = BoxTruncatedMvNormal(μₑ,Σₑ,a,b)
# # print(d)
# # moment(d,[0,1,1])
print(d)

# println("\n\n 1-d case:")

# μₑ = [3.5]
# Σₑ = 1*ones(1,1)
# a = [2.5]
# b = [4.5]
# d = BoxTruncatedMvNormal(μₑ,Σₑ,a,b)
# print(d)
# @show alpha(d)
# print(d)

# @time begin
#     a1 = -1.5; a2 = -1.4;
#     b1 = 1.7; b2 = 2.4;
#     μ₁ = 0.45; μ₂ = -0.2
#     σ₁ = 1.2; σ₂ = 0.8
#     ρ = 0.6
#     MMraw = [μ₁, μ₂]
#     SSraw = [σ₁^2 ρ*σ₁*σ₂ ;ρ*σ₁*σ₂ σ₂^2]
#     d = BoxTruncatedMvNormal(MMraw,SSraw,[a1,a2],[b1,b2])
# end

# print(d)