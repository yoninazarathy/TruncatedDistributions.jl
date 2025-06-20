using TruncatedDistributions
using Distributions
using HCubature
using PDMats
using LinearAlgebra

# using Pkg
# Pkg.activate(".")

# include("exampleDists.jl")

# μ = [0.5, -0.2]
# Σ = PDMat([1.5 0.0; 0.0 2.0])
# a1, b1 = -2, 3
# a2, b2 = -4, 1
# tnd1 = Truncated(Normal(μ[1],√Σ[1,1]), a1, b1)
# tnd2 = Truncated(Normal(μ[2],√Σ[2,2]), a2, b2)
# mnd1, vnd1 = mean(tnd1), var(tnd1)
# mnd2, vnd2 = mean(tnd2), var(tnd2)
# # @show (mnd1, vnd1)
# # @show (mnd2, vnd2)
# μ̂ = [mnd1, mnd2]
# Σ̂ = PDMat([vnd1 0.0; 0.0 vnd2])
# d = MvNormal(μ, Σ)
# r = BoxTruncationRegion([a1, a2], [b1, b2])
# dtrunc = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b; max_moment_levels = 4);
# μA =  mean(dtrunc)
# # @show μA
# grad = μ_gradient(dtrunc, μA, μ̂, Σ̂)


# # println("\nNow perturb μ")
# μ_off = μ + [0.2, 0.2]
# # @show norm(μ-μ_off)
# running_μ = copy(μ_off)

# if false #gradient descent just for μ
#     for i in 1:300
#         @show i   
#         global d = MvNormal(running_μ, Σ)
#         global r = BoxTruncationRegion([a1, a2], [b1, b2])
#         global dtrunc = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b; max_moment_levels = 4);
#         global μA = mean(dtrunc)
#         global grad = μ_gradient(dtrunc, μA, μ̂, Σ̂)
#         global running_μ = running_μ -0.2*grad
#         @show norm(running_μ - μ)
#     end
# end

# # println("\nCheck U gradient is zero at correct value")
# U = inv(Matrix(√Σ))
# μ_off = μ  + [-0.4, 0.3]
# @show U_gradient(dtrunc, μA, μ̂, Σ̂; U = U)
# @show U_gradient(dtrunc, μA, μ̂, Σ̂)

# println("\nNow perturb Σ")
# Σ_off = Σ + [0.15 0.04 ; 0.04 0.05]
# @show eigvals(Σ_off)
# @show norm(Σ - Σ_off)
# running_μ = copy(μ_off)
# running_U = cholesky(inv(Σ_off)).U
# running_Σ = inv(running_U)*inv(running_U)'
# # @show running_Σ, Σ_off

# errsΣ = []
# errsμ = []
# losses = []

# if false # complete gradient descent 
#     alpha = 0.01
#     println("\nGradient descent")
#     for i in 1:500
#         @show i
#         global d = MvNormal(running_μ, running_Σ)
#         global r = BoxTruncationRegion([a1, a2], [b1, b2])
#         global dtrunc = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b; max_moment_levels = 4);
#         global μA = mean(dtrunc)
#         global μ_grad = μ_gradient(dtrunc, μA, μ̂, Σ̂)
#         global U_grad = U_gradient(dtrunc, μA, μ̂, Σ̂)#QQQQ; U = running_U)
#         # if i == 70
#         #     global alpha = alpha/10
#         # end
#         global running_μ = running_μ - 5*alpha*μ_grad
#         global running_U = running_U - alpha*U_grad
#         # display(running_U)
#         global running_Σ = inv(running_U)*inv(running_U)'
#         display(running_Σ)
#         loss = norm(mean(dtrunc) - μ̂) + norm(cov(dtrunc) - Σ̂) 
#         push!(losses, loss)
#         # @show norm(running_μ - μ)
#         # @show norm(running_Σ - Σ)
#         # push!(errsμ, norm(running_μ - μ))
#         # push!(errsΣ, norm(running_Σ - Σ))
#     end
# end

# # hcub_m0 = hcubature((x)->pdf(dtrunc.untruncated,x),dtrunc.region.a, dtrunc.region.b)[1]
# # hcub_mx = hcubature((x)->pdf(dtrunc.untruncated,x).*x,dtrunc.region.a, dtrunc.region.b)[1]

# # U_gradient(dtrunc, μA, μ̂, Σ̂)

# Σ = [1 0.1 0.1;
#      0.1 2 0.1;
#      0.1 0.1 3]
# μ = [-1.5, 2.5, 1.3]
# d = MvNormal(μ, Σ)
# r = BoxTruncationRegion([-4, -2, -2], [3, 6, 4])
# dtrunc = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b; max_moment_levels = 4);
# @show mean(dtrunc)
# @show cov(dtrunc)
# @show tp(dtrunc)
# μ̂ = [-1.4, 2.6, 1.1]#mean(dtrunc)#[-1, 2.5, 3.3]
# Σ̂ = cov(dtrunc)
# PDMat([ 1.1  0.0  0.;
#         0.0  1.6  0.0;
#         0.00  0.0  1.7])
# #cov(dtrunc)

# Σ = PDMat([1 0.2;
#           0.2 2])
# μ = [2.5, 1.5]
# d = MvNormal(μ, Σ)
# r = BoxTruncationRegion([-2.8, -5], [4.25, 3.5])
# dtrunc = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b; max_moment_levels = 4);
# @show mean(dtrunc)
# @show cov(dtrunc)
# @show tp(dtrunc)
# μ̂ = [4, 2.]
# Σ̂ = PDMat([0.2 0.2;
#             0.2 0.8])
# # PDMat([0.4 0.15;
# #             0.15 1.3]);

# if false # complete gradient descent 
#     losses_μ = []
#     losses_Σ = []
#     dists = []
#     losses_total = []
#     running_μ = copy(μ)
#     running_Σ = copy(Σ)
#     running_U = cholesky(0.5*(inv(running_Σ) + inv(running_Σ)')).U
#     alpha = 0.005
#     println("\nGradient descent")
#     for i in 1:200
#         @show i
#         push!(dists, dtrunc)
#         global d = MvNormal(running_μ, 0.5*(running_Σ + running_Σ'))
#         global dtrunc = RecursiveMomentsBoxTruncatedMvNormal(d.μ, d.Σ, r.a, r.b; max_moment_levels = 4);
#         global μA = mean(dtrunc)
#         global μ_grad = μ_gradient(dtrunc, μA, μ̂, Σ̂)
#         global U_grad = U_gradient(dtrunc, μA, μ̂, Σ̂)
#         global running_μ = running_μ - 10*alpha*μ_grad
#         global running_U = running_U - alpha*U_grad
#         global running_Σ = inv(running_U)*inv(running_U)'
#         loss_μ = norm(mean(dtrunc) - μ̂)
#         loss_Σ =norm(cov(dtrunc) - Σ̂)
#         loss_total = loss_μ + loss_Σ
#         @show loss_total, (loss_μ, loss_Σ)
#         push!(losses_total, loss_total)
#         push!(losses_μ, loss_μ)
#         push!(losses_Σ, loss_Σ)
#     end
# end

# using Plots

# @info "Creating animation"
# @gif for dd in dists
#     x = range(r.a[1], r.b[1], length=300)
#     y = range(r.a[2], r.b[2], length=300)
#     X = repeat(x, 1, length(y))'
#     Y = repeat(y, 1, length(x))

#     # Compute the density at each grid point
#     Z = [log(pdf(dd.untruncated, [xi, yi])) for (xi, yi) in zip(vec(X), vec(Y))]
#     Z = reshape(Z, length(x), length(y))

#     # Plot the contours
#     contour(x, y, Z, 
#                 xlabel="x₁", 
#                 ylabel="x₂", 
#                 color=:viridis,
#                 xlim=(r.a[1], r.b[1]),
#                 ylim=(r.a[2], r.b[2]),
#                 aspect_ratio = 1, 
#                 legend=false,
#                 levels = 30)
# end

μ̂ = [4.9, -1.5]
Σ̂ = [7.5 0.3;
     0.3 0.4];
a = [-3.2, -3.5];
b = [4.5, 1.3];

dtrunc = loss_based_fit(μ̂, Σ̂, a, b)


# μ̂z = [0.0, 0]
# Σ̂z = [1  0.3/√(1.7*0.4);
#       0.3/√(1.7*0.4) 1]
# az = (a - μ̂) ./ diag(Σ̂)
# bz = (b - μ̂) ./ diag(Σ̂)

# dtruncz = loss_based_fit(μ̂z, Σ̂z, az, bz)
