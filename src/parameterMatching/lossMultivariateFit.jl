function loss_based_fit(μ̂::Vector{Float64}, Σ̂::Matrix{Float64}, a::Vector{Float64}, b::Vector{Float64}; 
                        μ_init::Vector{Float64} = μ̂, 
                        Σ_init::Matrix{Float64} = Σ̂,
                        α = 0.1,
                        min_grad_norm = 1e-4,
                        max_iter = Inf)
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

    μ_init0 = (μ_init- μ̂) ./ std_devs
    Σ_init0 = PDMat(Σ_init ./ (std_devs * std_devs'))


    losses_μ = []
    losses_Σ = []
    dists = []
    losses_total = []
    μ = μ_init0
    Σ = Σ_init0
    U = cholesky(0.5*(inv(Σ) + inv(Σ)')).U
    dtrunc = RecursiveMomentsBoxTruncatedMvNormal(μ, Σ, a0, b0; max_moment_levels = 4);

    
    grad_norm_sum = Inf
    i = 1
    while grad_norm_sum > min_grad_norm && i < max_iter
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
        loss_Σ = norm(cov(dtrunc) - Σ̂0)
        loss_total = loss_μ + loss_Σ
        grad_norm_sum = norm(μ_grad) + norm(U_grad)
        @show i, loss_total, grad_norm_sum#, (loss_μ, loss_Σ)
        push!(losses_total, loss_total)
        push!(losses_μ, loss_μ)
        push!(losses_Σ, loss_Σ)
        i+=1
    end
    Σ = PDMat(Σ .* (std_devs * std_devs')) #back to original coordinates
    μ = μ .*std_devs + μ̂
    return RecursiveMomentsBoxTruncatedMvNormal(μ, Σ, a, b; max_moment_levels = 4), 
                (losses_total = losses_total, 
                losses_μ = losses_μ,
                losses_Σ = losses_Σ,
                dists = dists)
end

function moment_loss(dist::RecursiveMomentsBoxTruncatedMvNormal, μ̂::AbstractVector{Float64}, Σ̂::AbstractMatrix{Float64})
    return 0.5*(norm(mean(dist) - μ̂)^2 + norm(cov(dist) - Σ̂)^2)
end

function vector_moment_loss(param_vec::Vector{Float64}, 
                            a,
                            b, 
                            μ̂::AbstractVector{Float64}, 
                            Σ̂::AbstractMatrix{Float64})
    μ, Σ = make_μ_Σ_from_param_vec(param_vec)
    dist = RecursiveMomentsBoxTruncatedMvNormal(μ, PDMat(Σ), a, b)
    return moment_loss(dist, μ̂, Σ̂)
end


function vector_gradient(   param_vec::Vector{Float64},
                            a,
                            b, 
                            μ̂::AbstractVector{Float64}, 
                            Σ̂::AbstractMatrix{Float64})
    μ, Σ = make_μ_Σ_from_param_vec(param_vec)
    dist = RecursiveMomentsBoxTruncatedMvNormal(μ, PDMat(Σ), a, b; max_moment_levels = 4)
    μA = mean(dist)
    μ_grad = μ_gradient(dist, μA, μ̂, PDMat(Σ̂))'
    U_grad = U_gradient(dist, μA, μ̂, PDMat(Σ̂))
    make_param_vec_from_μ_U(μ_grad, U_grad)
end

function make_μ_Σ_from_param_vec(param_vec)
    n = n_from_param_size(length(param_vec))
    μ = param_vec[1:n]
    U = [param_vec[3] param_vec[4]; 
         0.0          param_vec[5]] #do this for general n.
    Ui = inv(U)
    Σ = Ui*Ui'  #retrieve the covaranice from the upper triangular matrix
    return μ, Σ
end

function make_param_vec_from_μ_Σ(μ, Σ)
    F = cholesky(inv(Σ))
    U = F.U
    vcat(μ,U[1,1],U[1,2],U[2,2]) #first two coordinates are the mean and the last three coordinates are the factorized covariance
end

function make_param_vec_from_μ_U(μ, U)
    vcat(μ,U[1,1],U[1,2],U[2,2]) #first two coordinates are the mean and the last three coordinates are the factorized covariance
end

#n + n*(n+1)/2 = T
# 3n/2 + n^2/2  = T
# 3n + n^2 = 2T
# n^2 + 3n - 2T
# n = (-3 + sqrt(9 +8T))/2
function n_from_param_size(param_size::Integer)
    return Int((-3 + sqrt(9+8param_size))/2)
end

