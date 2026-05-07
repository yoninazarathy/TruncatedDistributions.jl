"""
    hcubature_inf(f, a, b; kwargs...)

Drop-in wrapper around `HCubature.hcubature` that accepts ±Inf entries in the
bound vectors. Each unbounded coordinate is mapped to a finite interval via
the standard substitutions (see https://github.com/stevengj/cubature#infinite-intervals):

* (-∞, +∞) :  x = t/(1−t²),       t ∈ (−1, 1),  dx/dt = (1+t²)/(1−t²)²
* ( a, +∞) :  x = a + t/(1−t),    t ∈ ( 0, 1),  dx/dt = 1/(1−t)²
* (-∞,  b) :  x = b − (1−t)/t,    t ∈ ( 0, 1),  dx/dt = 1/t²

The fully-finite case is dispatched straight to `HCubature.hcubature`.
"""
function hcubature_inf(f, a::AbstractVector, b::AbstractVector; kwargs...)
    n = length(a)
    @assert length(b) == n "bound vectors must have the same length"

    a_inf = [isinf(ai) && ai < 0 for ai in a]
    b_inf = [isinf(bi) && bi > 0 for bi in b]

    if !any(a_inf) && !any(b_inf)
        return hcubature(f, a, b; kwargs...)
    end

    a_t = Vector{Float64}(undef, n)
    b_t = Vector{Float64}(undef, n)
    for i in 1:n
        if a_inf[i] && b_inf[i]
            a_t[i], b_t[i] = -1.0, 1.0
        elseif a_inf[i]
            a_t[i], b_t[i] =  0.0, 1.0
        elseif b_inf[i]
            a_t[i], b_t[i] =  0.0, 1.0
        else
            a_t[i], b_t[i] = float(a[i]), float(b[i])
        end
    end

    function g(t)
        x   = Vector{Float64}(undef, n)
        jac = 1.0
        for i in 1:n
            ti = t[i]
            if a_inf[i] && b_inf[i]
                denom = 1 - ti^2
                x[i] = ti / denom
                jac *= (1 + ti^2) / denom^2
            elseif b_inf[i]              #  ( a, +∞ )
                denom = 1 - ti
                x[i] = a[i] + ti / denom
                jac *= 1 / denom^2
            elseif a_inf[i]              #  ( -∞, b )
                x[i] = b[i] - (1 - ti) / ti
                jac *= 1 / ti^2
            else
                x[i] = ti
            end
        end
        return f(x) * jac
    end

    return hcubature(g, a_t, b_t; kwargs...)
end
