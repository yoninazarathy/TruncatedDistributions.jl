"""
Compare four moment-matching methods on the canonical example set:

  1. full_gradient  – `loss_based_fit` (our explicit gradient)
  2. pair_gradient  – `pair_gradient_descent`  (block-of-2 coordinate descent on top of (1))
  3. optim_lbfgs    – `Optim.optimize` with LBFGS, finite-difference gradient
  4. prima_newuoa   – `PRIMA.newuoa`, derivative-free trust region

For each (example, method) we record (time, final moment loss, ok) and print a
summary at the end.

Run from the package root with:
    julia --project=. test/compare_methods.jl
"""

using TruncatedDistributions
using Distributions
using PDMats
using Printf

const METHODS = (
    full_gradient = correct_to_moments_with_full_gradient,
    pair_gradient = correct_to_moments_with_pair_gradient_descent,
    optim_lbfgs   = correct_to_moments_with_optim,
    prima_newuoa  = correct_to_moments_with_prima,
)

"""
Run `fit_fn(d, μ̂, Σ̂)`, returning (time_s, final_loss, ok). Errors are caught
so a single failure does not abort the sweep.
"""
function run_method(fit_fn, d, μ̂, Σ̂)
    try
        elapsed = @elapsed begin
            d_fit = fit_fn(d, μ̂, Σ̂)
        end
        return (time = elapsed, loss = moment_loss(d_fit, μ̂, Σ̂), ok = true)
    catch e
        @warn "method failed" exception = (e, catch_backtrace())
        return (time = NaN, loss = NaN, ok = false)
    end
end

"""
Build the moment-matching target by rounding the truncated moments of `d` to
`digits` decimal places. This guarantees a near-by feasible solution exists
while keeping the target slightly off from `(μ, Σ)`.
"""
function moment_target(d; digits::Int = 1)
    μ̂ = round.(mean(d); digits = digits)
    Σ̂ = round.(cov(d);  digits = digits)
    return μ̂, Σ̂
end

"""
Sweep all examples × all methods and return a Vector of NamedTuples, one per
(example, method) cell.
"""
function sweep(; digits::Int = 1)
    rows = NamedTuple[]
    for n in get_example_sizes()
        for i in 1:get_num_examples(n)
            d = dist_from_example(get_example(; n = n, index = i))
            μ̂, Σ̂ = moment_target(d; digits = digits)
            initial_loss = moment_loss(d, μ̂, Σ̂)
            @info "example" n = n index = i tp = round(tp(d), digits = 4) initial_loss = round(initial_loss, digits = 5)
            for (name, fn) in pairs(METHODS)
                r = run_method(fn, d, μ̂, Σ̂)
                push!(rows, (n = n, index = i, method = name, r...))
                @info "  $(rpad(string(name), 14))" time = round(r.time, digits = 3) loss = round(r.loss, digits = 5) ok = r.ok
            end
        end
    end
    return rows
end

function print_table(rows)
    println()
    @printf("%4s %5s %-14s %10s %12s %5s\n", "n", "idx", "method", "time(s)", "final_loss", "ok")
    println(repeat("-", 56))
    for r in rows
        @printf("%4d %5d %-14s %10.4f %12.6f %5s\n",
                r.n, r.index, string(r.method), r.time, r.loss, r.ok)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    rows = sweep()
    print_table(rows)
end
