using TruncatedDistributions
using Distributions
using HCubature
using PDMats
using LinearAlgebra
using Revise

μ̂ = [4.5, -1.0, 2.0]
Σ̂ = [0.3 0.2 -0.1;
     0.2 0.7 0.3;
     -0.1 0.3 1.2]
a = [3, -3.0, -4];
b = [6.0, 5, 6];
dtrunc_initial = RecursiveMomentsBoxTruncatedMvNormal(μ̂, PDMat(Σ̂),a,b)

@show tp(dtrunc_initial)
@show mean(dtrunc_initial)
@show cov(dtrunc_initial);

@info "pre-fitting"
mF, sF = zeros(3), zeros(3)
for i in 1:3
    mF[i], sF[i] = truncateDynamicFit(μ̂[i],sqrt(Σ̂[i,i]),(a[i],b[i]))
    @show i, mF[i], sF[i]^2
end

# mF2, sF2 =     truncateDynamicFit(μ̂[2],sqrt(Σ̂[2,2]),(a[2],b[2]))
# mF3, sF3 =     truncateDynamicFit(μ̂[3],sqrt(Σ̂[3,3]),(a[2],b[2]))
# @show (mF1,sF1^2)
# @show (mF2,sF2^2)
# μ0 = [mF1, mF2];
# Σ0 = [sF1^2 0 ; 0 sF2^2]


# using PRIMA
# @info "PRIMA"
# prima_result, prima_info = newuoa((v)->vector_moment_loss(v, a, b, μ̂, Σ̂), make_param_vec_from_μ_Σ(μ̂, Σ̂))
# @show prima_info
# μ_prima, Σ_prima = make_μ_Σ_from_param_vec(prima_result)
# dtrunc_prima = RecursiveMomentsBoxTruncatedMvNormal(μ_prima, PDMat(Σ_prima),a,b)
# @show mean(dtrunc_prima)
# @show cov(dtrunc_prima)
# @show tp(dtrunc_prima)

# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)
# @show mean(dtrunc)
# @show cov(dtrunc)