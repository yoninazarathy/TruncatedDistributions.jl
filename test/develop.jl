using TruncatedDistributions
using Distributions
using HCubature
using PDMats
using LinearAlgebra
using Revise

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




# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show cov(dtrunc)

# @info "Creating animation"
# @gif for (t,dd) in enumerate(logs.dists)
#      @show mean(dd)
#      @show cov(dd)
#     x = range(a[1], b[1], length=300)
#     y = range(a[2], b[2], length=300)
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
#                 xlim=(a[1], b[1]),
#                 ylim=(a[2], b[2]),
#                 aspect_ratio = 1, 
#                 legend=false,
#                 levels = 30, title = "t = $t, loss = $(round(logs.losses_total[t], digits = 5))")
# end every 1 fps=2




# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show cov(dtrunc)


# μ̂ = [4.5, -1.0, 2.0]
# Σ̂ = [0.8 0.0 0.0;
#      0.0 1.2 0.0;
#      0.0 0.0 0.5]
# a = [2, -4.0, -4];
# b = [6.5, 5, 5];

# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show cov(dtrunc)

# https://cran.r-project.org/web/packages/tmvtnorm/tmvtnorm.pdf and paper "Moments Calculation for the Doubly Truncated
# Multivariate Normal Density" (2021) B. G. Manjunath1 and Stefan Wilhelm2
# Example 1 in that paper
# μ̂ = [-0.1516343, -0.3881151]
# Σ̂ = [0.1630439 0.1613371;
#      0.1613371 0.6062505]
# a = [-1.0, -8]; #this -8 needs to be -inf but causes problems 
# b = [0.5, 1];

# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show cov(dtrunc)
# @show dtrunc.untruncated
# # interesting that we don't get the same original parameters in that paper, this is perhaps a situation of non-uniqueness.

# dtrunc_paper = RecursiveMomentsBoxTruncatedMvNormal([0.5, 0.5], PDMat([1.0 1.2 ; 1.2 2.0]), a, b);
# @show mean(dtrunc_paper)
# @show cov(dtrunc_paper)

# using Plots
# function plot_dd(dd,pp = plot(), color = :blue)
#      @show dd.untruncated
#     x = range(dd.region.a[1], dd.region.b[1], length=300)
#     y = range(max(dd.region.a[2], -3), dd.region.b[2], length=300)
#     X = repeat(x, 1, length(y))'
#     Y = repeat(y, 1, length(x))

#     # Compute the density at each grid point
#     Z = [log(pdf(dd.untruncated, [xi, yi])) for (xi, yi) in zip(vec(X), vec(Y))]
#     Z = reshape(Z, length(x), length(y))

#     # Plot the contours
#     return contour(pp,x, y, Z, 
#                     xlabel="x₁", 
#                     ylabel="x₂", 
#                     color=color,
#                     xlim=(a[1], b[1]),
#                     ylim=(max(a[2],-3), b[2]),
#                     aspect_ratio = 1, 
#                     legend=false,
#                     levels = 30,
#                     colorbar = false)
# end

# p = plot_dd(dtrunc)
# p = plot_dd(dtrunc_paper, p, :red)
# plot_dd(logs.dists[1], p, :black)

# # μ̂ = [-0.1516343, -0.3881151]
# Σ̂ = [0.1630439 0.0;
#      0.0 0.6062505]
# a = [-1.0, -8]; #this -8 needs to be -inf but causes problems 
# b = [0.5, 1];
# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show cov(dtrunc)

# μ̂ = [0.1]
# Σ̂ = fill(0.6^2, 1, 1)
# a = [-0.9]
# b = [1.35]
# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show sqrt(cov(dtrunc)[1][1])
# @show dtrunc.untruncated


# μ̂ = [-0.1516343, -0.3881151]
# Σ̂ = [0.1630439 0.1613371;
#      0.1613371 0.6062505]
# a = [-1.0, -8]; #this -8 needs to be -inf but causes problems 
# b = [0.5, 1];

# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show cov(dtrunc)
# @show dtrunc.untruncated

#first the uncorollated case
#better to do this with one-variable approach
# μ̂ = [-0.1516343, -0.3881151]
# Σ̂ = [0.1630439 0;
#      0 0.6062505]
# a = [-1.0, -5]; #this -3 needs to be -inf but causes problems 
# b = [0.5, 1];

# @info "First fit (independent)"
# ind_dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b, max_iter = 200)
# @show mean(ind_dtrunc)
# @show cov(ind_dtrunc)
# @show ind_dtrunc.untruncated.μ
# @show ind_dtrunc.untruncated.Σ

# @info "Second fit (independent)"
# Σ̂ = [0.1630439 0.1613371;
#      0.1613371 0.6062505]
# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b; μ_init = ind_dtrunc.untruncated.μ, Σ_init = ind_dtrunc.untruncated.Σ.mat)
# @show mean(dtrunc)
# @show cov(dtrunc)
# @show dtrunc.untruncated.μ
# @show dtrunc.untruncated.Σ

# using PDMats
μ̂ = [4.5, -1.0, 2.0]
Σ̂ = [0.8 0.1 -0.1;
     0.1 1.2 0.1;
     -0.1 0.1 0.5]
a = [2, -3.0, -4];
b = [6.5, 5, 2];
init_trunc = RecursiveMomentsBoxTruncatedMvNormal(μ̂, PDMat(Σ̂), a, b; max_moment_levels = 4);
@show mean(init_trunc)
@show tp(init_trunc)

# μ_g = μ_gradient(init_trunc, mean(init_trunc), μ̂, PDMat(Σ̂));
# U_g = U_gradient(init_trunc, mean(init_trunc), μ̂, PDMat(Σ̂));

# μ̂ = [0.3, -0.2]
# Σ̂ = [0.65 0.2; 
#      0.2 0.4]
# a = [-1.5, -1.4]
# b = [1.7, 2.4]
# @info "without pre-fitting"
# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show cov(dtrunc)
# @show tp(dtrunc)
# @show dtrunc.untruncated;


# @info "with pre-fitting"
# #Starting point for optimization comes from LiquetNazarathy paper (2015)
# mF1, sF1 =     truncateDynamicFit(μ̂[1],sqrt(Σ̂[1,1]),(a[1],b[1]))
# mF2, sF2 =     truncateDynamicFit(μ̂[2],sqrt(Σ̂[2,2]),(a[2],b[2]))
# @show (mF1,sF1^2)
# @show (mF2,sF2^2)
# μ0 = [mF1, mF2];
# Σ0 = [sF1^2 0 ; 0 sF2^2]

# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b;μ_init = μ0, Σ_init = Σ0)
# @show mean(dtrunc)
# @show cov(dtrunc)
# @show tp(dtrunc)
# @show dtrunc.untruncated;

# μ̂ = [4.5, -1.0]

# Σ̂ = [0.8 0.3
#      0.3 0.2];
# a = [3.5, -2.0];
# b = [6.5, 1];
# d_init = RecursiveMomentsBoxTruncatedMvNormal(μ̂, PDMat(Σ̂),a,b)
# μA = mean(d_init)
# @show tp(d_init)

# using PRIMA
# @info "PRIMA"
# @time prima_result, prima_info = newuoa((v)->vector_moment_loss(v, a, b, μ̂, Σ̂), make_param_vec_from_μ_Σ(μ̂, Σ̂))
# @show prima_info
# μ_prima, Σ_prima = make_μ_Σ_from_param_vec(prima_result)
# dtrunc_prima = RecursiveMomentsBoxTruncatedMvNormal(μ_prima, PDMat(Σ_prima),a,b)
# @show mean(dtrunc_prima)
# @show cov(dtrunc_prima)
# @show tp(dtrunc_prima)

using Optim
@info "no gradient supplied - optim"
# optim_result = optimize((v)->vector_moment_loss(v, a, b, μ̂, Σ̂),  #function
#                         make_param_vec_from_μ_Σ(μ̂, Σ̂), #initial value
#                         LBFGS(),
#                         Optim.Options(show_trace=true, time_limit=10))# DOESN'T WORK:,autodiff = :forward)
# μ_optim, Σ_optim = make_μ_Σ_from_param_vec(optim_result.minimizer)
# dtrunc_optim = RecursiveMomentsBoxTruncatedMvNormal(μ_optim, PDMat(Σ_optim),a,b)
# @show mean(dtrunc_optim)
# @show cov(dtrunc_optim)

# using Optim
# # @info "gradient Optim - no standardization" 
# init_v = make_param_vec_from_μ_Σ(μ̂, Σ̂)
# for i in 1:3
#   @info "Iteration $i"
#   global μA, init_v
#   @show init_v
#   optim_grad_result = optimize((v)->approximate_vector_moment_loss(v, a, b, μA, μ̂, Σ̂),  #function
#                                   (v)->vector_gradient(v, a, b, μA, μ̂, Σ̂), #gradient
#                                   init_v, #initial value
#                                   LBFGS(), #doesn't work 
#                                   # Adam(alpha=0.01),#, beta_mean=0.0, beta_var=0.0),
#                                   Optim.Options(show_trace=true, iterations = 10),#g_calls_limit = 1), #time_limit =50),#, f_abstol = 1e-6),
#                                 inplace = false)
#   μ_optim_grad, Σ_optim_grad = make_μ_Σ_from_param_vec(optim_grad_result.minimizer)
#   dtrunc_optim_grad = RecursiveMomentsBoxTruncatedMvNormal(μ_optim_grad, PDMat(Σ_optim_grad),a,b)
#   @show moment_loss(dtrunc_optim_grad, μ̂, Σ̂)
#   @show approximate_moment_loss(dtrunc_optim_grad, μA, μ̂, Σ̂)
#   @show μA
#   μA = mean(dtrunc_optim_grad)
#   @show μA
#   init_v = make_param_vec_from_μ_Σ(μ_optim_grad, Σ_optim_grad)
#   @show approximate_vector_moment_loss(init_v, a, b, μA, μ̂, Σ̂)
# end
# @show mean(dtrunc_optim_grad)
# @show cov(dtrunc_optim_grad)


# @info "standardization"
# std_devs = sqrt.(diag(Σ̂))
# μ̂0 = zeros(2);
# Σ̂0 =  PDMat(Σ̂ ./ (std_devs * std_devs'))
# a0 = (a - μ̂) ./ std_devs
# b0 = (b - μ̂) ./ std_devs
# prima_result0, prima_info0 = newuoa((v)->vector_moment_loss(v, a0, b0, μ̂0, Σ̂0), make_param_vec_from_μ_Σ(μ̂0, Σ̂0))
# # FAILS (or maybe not - check-above)
# optim_result0 = optimize((v)->vector_moment_loss(v, a0, b0, μ̂0, Σ̂0), 
#                         make_param_vec_from_μ_Σ(μ̂0, Σ̂0),
#                         LBFGS(),        
#                         Optim.Options(show_trace=true,time_limit=10))#f_abstol = 1e-5))
# μ_optim0, Σ_optim0 = make_μ_Σ_from_param_vec(optim_result0.minimizer)
# dtrunc_optim0 = RecursiveMomentsBoxTruncatedMvNormal(μ_optim0, PDMat(Σ_optim0),a0,b0)
# @show mean(dtrunc_optim0)
# @show cov(dtrunc_optim0)
# μ_optim = μ_optim0 .*std_devs + μ̂
# Σ_optim = PDMat(Σ_optim0 .* (std_devs * std_devs')) #back to original coordinates
# dtrunc_optim = RecursiveMomentsBoxTruncatedMvNormal(μ_optim, PDMat(Σ_optim),a,b)
# @show mean(dtrunc_optim)
# @show cov(dtrunc_optim)
# @show tp(dtrunc_optim0)
# @show tp(dtrunc_optim)


# std_devs = sqrt.(diag(Σ̂))
# μ̂0 = zeros(2);
# Σ̂0 =  PDMat(Σ̂ ./ (std_devs * std_devs'))
# a0 = (a - μ̂) ./ std_devs
# b0 = (b - μ̂) ./ std_devs
# optim_result0 = optimize((v)->vector_moment_loss(v, a0, b0, μ̂0, Σ̂0), 
#                          (v)->vector_gradient(v, a, b, μ̂0, Σ̂0),
#                                 make_param_vec_from_μ_Σ(μ̂0, Σ̂0),
#                               #   LBFGS(),
#                                 Adam(alpha=0.01),# beta_mean=0.0, beta_var=0.0),
#                                 # GradientDescent(alpha=0.1),
#                                 Optim.Options(show_trace=true, f_abstol = 1e-6),
#                                 inplace = false)
# μ_optim0, Σ_optim0 = make_μ_Σ_from_param_vec(optim_result0.minimizer)
# dtrunc_optim0 = RecursiveMomentsBoxTruncatedMvNormal(μ_optim0, PDMat(Σ_optim0),a0,b0)
# @show mean(dtrunc_optim0)
# @show cov(dtrunc_optim0)
# μ_optim = μ_optim0 .*std_devs + μ̂
# Σ_optim = PDMat(Σ_optim0 .* (std_devs * std_devs')) #back to original coordinates
# dtrunc_optim = RecursiveMomentsBoxTruncatedMvNormal(μ_optim, PDMat(Σ_optim),a,b)
# @show mean(dtrunc_optim)
# @show cov(dtrunc_optim)
# @show tp(dtrunc_optim0)
# @show tp(dtrunc_optim)


# Just playing with LBFGS
# @info "no gradient"
# optim_result = optimize((x)->(x-[1,1])'*(x-[1,1]),[10.0,10],LBFGS(),Optim.Options(show_trace=true)) # no graident
# @show optim_result.minimizer

# @info "yes gradient"
# optim_result = optimize((x)->(x-[1,1])'*(x-[1,1]),(x)->2(x .- 1), [10.,10], LBFGS(),Optim.Options(show_trace=true), inplace = false)
# @show optim_result.minimizer

# @info "yes gradient - but slightly wrong"
# optim_result = optimize((x)->(x-[1,1])'*(x-[1,1]),(x)->2.2(x .- 1), [10.,10], LBFGS(),Optim.Options(show_trace=true), inplace = false)
# @show optim_result.minimizer