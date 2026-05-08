"""
Experiment (Phase 2a): hand LBFGS the *matched* (loss, gradient) pair.

Hypothesis: previous LBFGS attempts failed because the loss function fed in
was the true `moment_loss` while the gradient was `vector_gradient`, which is
the gradient of the *surrogate* `approximate_moment_loss` (with μA frozen).
Mismatched f and ∇f kill LBFGS's line search.

This script holds μA = μ̂ fixed and gives LBFGS the matched pair:
    f(p)    = approximate_vector_moment_loss(p; a, b, μA=μ̂, μ̂, Σ̂)
    ∇f(p)   = vector_gradient(p; a, b, μA=μ̂, μ̂, Σ̂)

Run from the package root with:
    julia --project=. test/experiment_lbfgs_our_gradient.jl
"""

using TruncatedDistributions
using Distributions, PDMats, LinearAlgebra
using Optim
using Printf

function lbfgs_with_our_gradient(d, μ̂::Vector{Float64}, Σ̂::Matrix{Float64};
                                  μA_fixed = nothing, max_iter = 200)
    a, b = d.region.a, d.region.b
    μA = μA_fixed === nothing ? copy(μ̂) : copy(μA_fixed)

    f(p)     =  approximate_vector_moment_loss(p, a, b, μA, μ̂, Σ̂)
    g!(g, p) = (g .= vector_gradient(p, a, b, μA, μ̂, Σ̂))

    p0 = make_param_vec_from_μ_Σ(μ̂, Σ̂)

    res = optimize(f, g!, p0, LBFGS(),
                   Optim.Options(iterations = max_iter,
                                 g_tol      = 1e-6,
                                 show_trace = false))
    μ_fit, Σ_fit = make_μ_Σ_from_param_vec(res.minimizer)
    d_fit = RecursiveMomentsBoxTruncatedMvNormal(μ_fit, PDMat(Σ_fit), a, b)
    return d_fit, res
end

function run_one(label::String, ne)
    d  = dist_from_example(ne)
    μ̂  = collect(ne.μ̂)
    Σ̂  = Matrix(ne.Σ̂)

    @info "=== $label ===" μ=ne.μ Σ=ne.Σ a=ne.a b=ne.b μ̂=μ̂ Σ̂=Σ̂

    p0     = make_param_vec_from_μ_Σ(μ̂, Σ̂)
    a, b   = d.region.a, d.region.b
    μA     = copy(μ̂)
    f0     = approximate_vector_moment_loss(p0, a, b, μA, μ̂, Σ̂)
    true0  = moment_loss(RecursiveMomentsBoxTruncatedMvNormal(μ̂, PDMat(Σ̂), a, b), μ̂, Σ̂)
    @info "initial" surrogate_loss = f0 true_loss = true0

    elapsed = @elapsed begin
        d_fit, res = lbfgs_with_our_gradient(d, μ̂, Σ̂)
    end
    pf       = make_param_vec_from_μ_Σ(d_fit.untruncated.μ, Matrix(d_fit.untruncated.Σ))
    fend     = approximate_vector_moment_loss(pf, a, b, μA, μ̂, Σ̂)
    trueend  = moment_loss(d_fit, μ̂, Σ̂)

    @info "result" iters = res.iterations f_calls = res.f_calls g_calls = res.g_calls converged = Optim.converged(res) elapsed = round(elapsed, digits=3)
    @info "loss"   surrogate_final = fend true_final = trueend
    @info "iterate" μ_fit = d_fit.untruncated.μ Σ_fit = Matrix(d_fit.untruncated.Σ)
    println()
    return (label = label, time = elapsed,
            surrogate_init = f0, surrogate_final = fend,
            true_init = true0, true_final = trueend,
            iters = res.iterations, converged = Optim.converged(res))
end

function main()
    rows = []
    push!(rows, run_one("Example 1 (MW2021, semi-infinite)", get_example(n=2, index=3)))
    push!(rows, run_one("Example 2 (fully bounded box)",     get_example(n=2, index=4)))

    println()
    @printf("%-44s %8s %12s %12s %12s %12s %5s\n",
            "case", "time(s)", "surr_init", "surr_final", "true_init", "true_final", "conv")
    println(repeat("-", 110))
    for r in rows
        @printf("%-44s %8.3f %12.4e %12.4e %12.4e %12.4e %5s\n",
                r.label, r.time, r.surrogate_init, r.surrogate_final,
                r.true_init, r.true_final, r.converged)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
