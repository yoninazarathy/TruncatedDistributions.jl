function correct_to_moments_with_prima(     d::RecursiveMomentsBoxTruncatedMvNormal, 
                                            μ̂::AbstractVector{Float64},
                                            Σ̂::AbstractMatrix{Float64})
    prima_result, prima_info = newuoa((v)->vector_moment_loss(v, d.region.a, d.region.b, μ̂, Σ̂), make_param_vec_from_μ_Σ(μ̂, Σ̂))
    # @show prima_info
    μ_prima, Σ_prima = make_μ_Σ_from_param_vec(prima_result)
    return RecursiveMomentsBoxTruncatedMvNormal(μ_prima, PDMat(Σ_prima),d.region.a, d.region.b)
end

function correct_to_moments_with_optim(     d::RecursiveMomentsBoxTruncatedMvNormal, 
                                            μ̂::AbstractVector{Float64},
                                            Σ̂::AbstractMatrix{Float64})
    optim_result = optimize((v)->vector_moment_loss(v, d.region.a, d.region.b, μ̂, Σ̂),  #function
                            make_param_vec_from_μ_Σ(μ̂, Σ̂), #initial value
                            LBFGS(),
                            Optim.Options(show_trace=false))
    μ_optim, Σ_optim = make_μ_Σ_from_param_vec(optim_result.minimizer)
    return RecursiveMomentsBoxTruncatedMvNormal(μ_optim, PDMat(Σ_optim),d.region.a, d.region.b)
end