using TruncatedDistributions
using Distributions
using HCubature
using PDMats
using LinearAlgebra
using Revise
using Plots

# μ̂ = [4.5, -1.0, 2.0]
# Σ̂ = [0.3 0.2 -0.1;
#      0.2 0.7 0.3;
#      -0.1 0.3 1.2]
# a = [0, -3.0, -6];
# b = [6.0, 5, 6];
# dtrunc_initial = RecursiveMomentsBoxTruncatedMvNormal(μ̂, PDMat(Σ̂),a,b)

# @show tp(dtrunc_initial)
# @show mean(dtrunc_initial)
# @show cov(dtrunc_initial);

# dtrunc, logs = loss_based_fit(μ̂, Σ̂, a, b)

# μ_running = copy(μ̂) 
# Σ_running = copy(Σ̂)

# Subproblem just on first two coordinates 
# current_set = [1,2]
# μ̂1 = μ̂[current_set]
# Σ̂1 = Σ̂[current_set,current_set]
# a1 = a[current_set]
# b1 = b[current_set];
#fit a 2x2 problem on the current set
# dtrunc, logs = loss_based_fit(μ̂1, Σ̂1, a1, b1; 
#                                 μ_init = μ_running[current_set],
#                                 Σ_init = Σ_running[current_set, current_set],
#                                 max_iter=10);
# dtrunc.untruncated

# μ_running[current_set] = dtrunc.untruncated.μ
# Σ_running[current_set, current_set] = dtrunc.untruncated.Σ

# Second iteration Now on coordinates 1 and 3 
# current_set = [1,3]
# μ̂2 = μ̂[current_set]
# Σ̂2 = Σ̂[current_set,current_set]
# a2 = a[current_set]
# b2 = b[current_set];

# dtrunc, logs = loss_based_fit(μ̂2, Σ̂2, a2, b2;
#                                 μ_init = μ_running[current_set],
#                                 Σ_init = Σ_running[current_set, current_set],
#                                 max_iter=10);
# dtrunc.untruncated

μ̂ = [4.5, -1.0, 2.0]
Σ̂ = [0.3 0.2 -0.1;
     0.2 0.7 0.3;
     -0.1 0.3 1.2]
a = [1, -3.0, -6];
b = [5.5, 4.5, 6];
dtrunc_not_adjusted = RecursiveMomentsBoxTruncatedMvNormal(μ̂, PDMat(Σ̂),a,b);
@show tp(dtrunc_not_adjusted)
@show mean(dtrunc_not_adjusted)
@show cov(dtrunc_not_adjusted);

μ_running = copy(μ̂) 
Σ_running = copy(Σ̂)

total_loss = []

possible_sets = [[1,2], [1,3], [2,3]]

function find_set_with_worst_loss()
    set_losses = zeros(length(possible_sets))
    for (i,s) in enumerate(possible_sets)
        μ_s = μ_running[s]
        Σ_s = Σ_running[s,s]
        μ̂_s = μ̂[s]
        Σ̂_s = Σ̂[s,s]
        a_s = a[s]
        b_s = b[s];
        dtrunc = RecursiveMomentsBoxTruncatedMvNormal(μ_s, PDMat(Σ_s),a_s,b_s)
        all_coordinates_loss = moment_loss(dtrunc, μ̂_s, Σ̂_s)
        set_losses[i] = all_coordinates_loss
    end
    @show set_losses
    return possible_sets[argmax(set_losses)]
end

for i in 1:150
    current_set = find_set_with_worst_loss()
    @info "Starting iteration $i on set $current_set"

    #create the 2x2 problem
    μ̂_current = μ̂[current_set]
    Σ̂_current = Σ̂[current_set,current_set]
    a_current = a[current_set]
    b_current = b[current_set];

    #do a few gradient descent steps on the 2x2 problem
    dtrunc, logs = loss_based_fit(μ̂_current, Σ̂_current, a_current, b_current; 
                                    μ_init = μ_running[current_set],
                                    Σ_init = Σ_running[current_set, current_set],
                                    max_iter=3, α = 0.01);
    
    #update the over all problem
    μ_running[current_set] = dtrunc.untruncated.μ
    Σ_running[current_set, current_set] = dtrunc.untruncated.Σ

    display(μ_running)
    display(Σ_running)

    dtrunc_all_coords = RecursiveMomentsBoxTruncatedMvNormal(μ_running, PDMat(Σ_running),a,b)
    all_coordinates_loss = moment_loss(dtrunc_all_coords, μ̂, Σ̂)
    @show all_coordinates_loss
    push!(total_loss, all_coordinates_loss)
end

# @info "finished"
# @show mean(dtrunc_all_coords)
# @show cov(dtrunc_result)

plot(total_loss)
