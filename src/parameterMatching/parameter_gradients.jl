function μ_gradient(    d::RecursiveMomentsBoxTruncatedMvNormal, 
                        μA::Vector{Float64}, 
                        μ̂::Vector{Float64}, 
                        Σ̂::PDMat)
    μ = d.untruncated.μ
    Σ = d.untruncated.Σ
    n = length(d)
    m(inds) = raw_moment_from_indices(d,inds)
    m0 = m(Int[])

    ii(k) = ((k-1) % n)+1
    jj(k) = floor(Int,(k-1)/n)+1

    I1 = [m([i]) - m0*μ[i] for i in 1:n]
    I2 = [m([i,j]) - m([i])*μ[j] - m([j])*μ̂[i] + m0*μ̂[i]*μ[j] for i in 1:n, j in 1:n]
    I3 = [m([ii(k), jj(k)]) - m([ii(k)])*μA[jj(k)] - m([jj(k)])*μA[ii(k)] + m0*(μA[ii(k)])*μA[jj(k)] - Σ̂[ii(k),jj(k)] for k in 1:n^2]
    I4 = [m([ii(k),jj(k),l]) - m([ii(k),l])*μA[jj(k)] - m([jj(k),l])*μA[ii(k)] - m([ii(k),jj(k)])*μ[l] + m([l])*μA[ii(k)]*μA[jj(k)] + m([jj(k)])*μA[ii(k)]*μ[l] - m0*(μA[ii(k)]*μA[jj(k)] - Σ̂[ii(k),jj(k)]) for k in 1:n^2, l in 1:n]

    return (I1'*I2 + I3'*I4)*inv(Σ)  #QQQQ check inv... 
end

function U_gradient(    d::RecursiveMomentsBoxTruncatedMvNormal, 
                        μA::Vector{Float64}, 
                        μ̂::Vector{Float64}, 
                        Σ̂::PDMat)
    μ = d.untruncated.μ
    Σ = d.untruncated.Σ
    n = length(d)

    m(inds) = raw_moment_from_indices(d,inds)
    m0 = m(Int[])

    I1(i) = fill(m([i]) - m0*μ[i],n,n)
    I2(i) = [k < l ? 0 : m([i,k,l]) + m([i,k])*μ[l] + m([k,l])*μ[i] + m([i,l])*μ[k] + m([i])*(Σ[k,l] + μ[k]*μ[l]) + m([l])*μ[i]*μ[k] + m([k])*μ[i]*μ[l] - m0*μ[i]*(Σ[k,l] - μ[k]*μ[l]) for k in 1:n, l in 1:n]
    I3(i) = 


    sum([I1(i) .* I2(i) for i in 1:n]) #QQQQ missing U multiplications
end