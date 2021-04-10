cd(@__DIR__)
using Pkg; Pkg.activate("..")
using Revise
# Pkg.precompile()

using TruncatedDistributions
using Test
using Distributions
using HCubature

# a = [-20.,-20,-20,-20,-20]
# b = [20.,20.,20,20,20]

# r = BoxTruncationRegion(a,b)
# intruncationregion(r,[0.0,18,20,0,0])
# d = MvN

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

# μₑ = [2.5, 3.5]
# Σₑ = [2.0 -0.5;
#        -0.5 5.0]
# a = [-2.3,-20.]
# b = [4.,12.3]


# untruncated = MvNormal(μₑ,Σₑ)
# region = BoxTruncationRegion(a,b)
# d = TruncatedMvDistribution(untruncated,region)

# x = [10.,2]
# @show insupport(d,x)
# @show rand(d)
# @show pdf(d,x)
# @show mean(d)
# @show cov(d)
# # compute_tp(d)
# # compute_mean(d)
# # compute_cov(d)
# @show d.state.tp
# @show d.state.μ
# @show d.state.Σ

# # moment(7)

# td = TruncatedNormal(0.25,1.0,-1.,1.)
# @show moment(td,2)
# @show mean(td)^2 + var(td)
# # @show moments(td,22)
# td

# @show moment(td,10)

# d = BoxTruncatedMvNormalRecursiveMoments(μₑ,Σₑ,a,b)
# # # print(d)
# # # moment(d,[0,1,1])
# print(d)

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
#     d = BoxTruncatedMvNormalRecursiveMoments(MMraw,SSraw,[a1,a2],[b1,b2],4)
# end

# untruncated = MvNormal(MMraw,SSraw)
# region = BoxTruncationRegion([a1,a2],[b1,b2])
# d2 = TruncatedMvDistribution(untruncated,region)

# μₑ = [3.5,2,3.5,3.5]
# Σₑ = [ 7. 1 0 1 ; 
#        1 3.3 2 0 ;
#        0 2 3.4 0 ;
#        1 0 0 4 ]
# a = [-20. ,-20, -20 ,-20]
# b = [20. ,20. ,20, 20]

# dist1 = BoxTruncatedMvNormalRecursiveMoments(μₑ,Σₑ,a,b,2)

# @show cov(dist1)

# untruncated = MvNormal(μₑ,Σₑ)
# region = BoxTruncationRegion(a,b)
# dist2 = TruncatedMvDistribution(untruncated,region)
# @show cov(dist2)

μₑ = [3.5,2,3.5]
Σₑ = [ 7. 1 0 ; 
       1 3.3 2  ;
       0 2 3.8 ]
a = [-2. ,-1 ,2]
b = [4. ,4. , 5]

# dist1 = BoxTruncatedMvNormalRecursiveMoments(μₑ,Σₑ,a,b,2)

# @show cov(dist1)

# untruncated = MvNormal(μₑ,Σₑ)
# region = BoxTruncationRegion(a,b)
# dist2 = TruncatedMvDistribution{MvNormal,BoxTruncationRegion,TruncatedMvDistributionSecondOrderState}(untruncated,region)#,state)
# @show cov(dist2)
d = BasicBoxTruncatedMvNormal(μₑ,Σₑ,a,b)
d2 = RecursiveMomentsBoxTruncatedMvNormal(μₑ,Σₑ,a,b)
compute_moments(d2.state)

# mean(d2)
# abstract type Foo end

# struct ConcreteFoo <: Foo
#         x::Int
# end

# struct Bar{T <: Real, F <: Foo}
#         u::Vector{T}
#         f::F
#         # function Bar(v::Vector{Float64},f::F) where F <: Foo
#         #         new{F}(v,f)
#         # end
# end

# function Bar{T,F}(v::Vector{T}) where {T <: Real, F <: Foo}
#         Bar(v,F(length(v)))
# end

# Bar{Float64,ConcreteFoo}([2.3,5.2])


# struct Foo
#         x::Float64
#         Foo{Int}() = new{typeof{Foot{Int}}}(3.14)
# end

# object = Foo{5}()