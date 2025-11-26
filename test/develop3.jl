using TruncatedDistributions
using Distributions
using PDMats

for n in [4]#get_example_sizes()
    for i in 1:get_num_examples(n)
        @info (n, i)
        d = dist_from_example(get_example(;n=n,index=i))
        μ̂ = round.(mean(d), digits = 1)
        Σ̂ = round.(cov(d), digits = 1)
        @show initial_moment_loss = moment_loss(d, μ̂, Σ̂)

        @info "Pair gradient descent based fit"
        pair_gradient_based_not_fail = true
        try
            time_pair_gradient_based = @elapsed begin
                d_pair_gradient_based = correct_to_moments_with_pair_gradient_descent(d, μ̂, Matrix(Σ̂))
            end
            pair_gradient_based_moment_loss = moment_loss(d_pair_gradient_based, μ̂, Σ̂)
            @show time_pair_gradient_based, pair_gradient_based_moment_loss, pair_gradient_based_not_fail
        catch e
            println("Caught an error", e)
            pair_gradient_based_not_fail = false
        end

        # @info "Gradient based fit"
        # gradient_based_not_fail = true
        # try
        #     time_gradient_based = @elapsed begin
        #         d_gradient_based = loss_based_fit(μ̂, Matrix(Σ̂), d.region.a, d.region.b)
        #     end
        #     gradient_based_moment_loss = moment_loss(d_gradient_based, μ̂, Σ̂)
        #     @show time_gradient_based, gradient_based_moment_loss, gradient_based_not_fail
        # catch e
        #     println("Caught an error")#, e)
        #     gradient_based_not_fail = false
        # end

        # @info "Optim"
        # optim_not_fail = true
        # try
        #     time_optim = @elapsed begin
        #         d_optim = correct_to_moments_with_optim(d, μ̂, PDMat(Σ̂))
        #     end
        #     optim_moment_loss = moment_loss(d_optim, μ̂, Σ̂)
        #     @show time_optim, optim_moment_loss, optim_not_fail
        # catch e
        #     println("Caught an error")#, e)
        #     optim_not_fail = false

        # end

        # @info "PRIMA"
        # prima_not_fail = true
        # try
        #     time_prima = @elapsed begin
        #         d_prima = correct_to_moments_with_prima(d, μ̂, PDMat(Σ̂))
        #     end
        #     prima_moment_loss = moment_loss(d_prima, μ̂, Σ̂)
        #     @show time_prima, prima_moment_loss, prima_not_fail
        # catch e
        #     println("Caught an error")#, e)
        #     prima_not_fail = false
        # end

    end
end