
function raw_moment(d::RecursiveMomentsBoxTruncatedMvNormal, k::Vector{Int}) 
    raw_moment(d.state, k)
end

function raw_moment_dict(d::RecursiveMomentsBoxTruncatedMvNormal) 
    raw_moment_dict(d.state)
end

function raw_moment_from_indices(d::RecursiveMomentsBoxTruncatedMvNormal, indices::Vector{Int}) 
    kappa = zeros(Int, length(d))
    for i in indices
        kappa[i] += 1
    end
    raw_moment(d, kappa)
end