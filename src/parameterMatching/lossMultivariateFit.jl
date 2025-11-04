function loss_based_fit(μ̂::Vector{Float64}, Σ̂::Matrix{Float64}, a::Vector{Float64}, b::Vector{Float64})
    n = length(μ̂)
    n1, m1 = size(Σ̂)
    (n != n1 || n != m1) && error("Mismatch of dimensions")  
    
    std_devs = sqrt.(diag(Σ̂))
    μ̂0 = zeros(n);
    Σ̂0 =  PDMat(Σ̂ ./ (std_devs * std_devs'))
    a0 = (a - μ̂) ./ std_devs
    b0 = (b - μ̂) ./ std_devs
    !all(a0 .< 0.0) && error("a needs to be less than μ̂")
    !all(b0 .> 0.0) && error("b needs to be less than μ̂")

    #QQQQ check if sigma_i is too high for a given [a,b] then flag it is not possible - error 

    losses_μ = []
    losses_Σ = []
    dists = []
    losses_total = []
    μ = copy(μ̂0)
    Σ = PDMat(copy(Σ̂0))
    U = cholesky(0.5*(inv(Σ) + inv(Σ)')).U
    dtrunc = RecursiveMomentsBoxTruncatedMvNormal(μ, Σ, a0, b0; max_moment_levels = 4);

    α = 0.1
    @showprogress for _ in 1:5000
        push!(dists,  RecursiveMomentsBoxTruncatedMvNormal(dtrunc.untruncated.μ .*std_devs + μ̂, PDMat(dtrunc.untruncated.Σ .* (std_devs * std_devs')), a, b; max_moment_levels = 2))
        dtrunc = RecursiveMomentsBoxTruncatedMvNormal(μ, Σ, a0, b0; max_moment_levels = 4);
        μA = mean(dtrunc)
        μ_grad = μ_gradient(dtrunc, μA, μ̂0, Σ̂0)'
        U_grad = U_gradient(dtrunc, μA, μ̂0, Σ̂0) #QQQQ pass in U
        μ = μ - 10*α*μ_grad #gradient based update
        U = U - α*U_grad #gradient based update
        Ui = inv(U)
        Σ = PDMat(Ui*Ui')
        loss_μ = norm(mean(dtrunc) - μ̂0)
        loss_Σ =norm(cov(dtrunc) - Σ̂0)
        loss_total = loss_μ + loss_Σ
        @show loss_total, (loss_μ, loss_Σ)
        push!(losses_total, loss_total)
        push!(losses_μ, loss_μ)
        push!(losses_Σ, loss_Σ)
    end
    Σ = PDMat(Σ .* (std_devs * std_devs')) #back to original coordinates
    μ = μ .*std_devs + μ̂
    return RecursiveMomentsBoxTruncatedMvNormal(μ, Σ, a, b; max_moment_levels = 4), 
                (losses_total = losses_total, 
                losses_μ = losses_μ,
                losses_Σ = losses_Σ,
                dists = dists)
end