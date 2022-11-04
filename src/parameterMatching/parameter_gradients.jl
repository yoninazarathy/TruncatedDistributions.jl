function μ_gradient(    d::RecursiveMomentsBoxTruncatedMvNormal, 
                        μA::Vector{Float64}, 
                        μ̂::Vector{Float64}, 
                        Σ̂::PDMat)
    μ = d.untruncated.μ
    Σ = d.untruncated.Σ
    n = length(d)
    m(inds) = raw_moment_from_indices(d, inds)
    m0 = m(Int[])

    ii(k) = ((k-1) % n) + 1
    jj(k) = floor(Int,(k-1)/n) + 1

    I1 = [m([i]) - m0*μ̂[i] for i in 1:n]
    I2 = [m([i,j]) - m([i])*μ[j] - m([j])*μ̂[i] + m0*μ̂[i]*μ[j] for i in 1:n, j in 1:n]
    I3 = [m([ii(k), jj(k)]) - m([ii(k)])*μA[jj(k)] - m([jj(k)])*μA[ii(k)] + m0*(μA[ii(k)])*μA[jj(k)] - Σ̂[ii(k),jj(k)] for k in 1:n^2]
    I4 = [m([ii(k),jj(k),l]) - m([ii(k),l])*μA[jj(k)] - m([jj(k),l])*μA[ii(k)] - m([ii(k),jj(k)])*μ[l] + m([l])*μA[ii(k)]*μA[jj(k)] + m([jj(k)])*μA[ii(k)]*μ[l] - m0*(μA[ii(k)]*μA[jj(k)] - Σ̂[ii(k),jj(k)]) for k in 1:n^2, l in 1:n]

    return (I1'*I2 + I3'*I4)*inv(Σ) 
end

function U_gradient(    d::RecursiveMomentsBoxTruncatedMvNormal, 
                        μA::Vector{Float64}, 
                        μ̂::Vector{Float64}, 
                        Σ̂::PDMat;
                        U::Matrix{Float64} #QQQQ recomendation...)
    μ = d.untruncated.μ
    Σ = d.untruncated.Σ
    U =  .....QQQQ compute U if it is not given....
    n = length(d)

    m(inds) = raw_moment_from_indices(d, inds)
    m0 = m(Int[])

    I1(i) = fill(m([i]) - m0*μ̂[i], n, n)
    I2(i) = [k < l ? 0 : m([i,k,l]) + m([i,k])*μ[l] + m([k,l])*μ̂[i] + m([i,l])*μ[k] + m([i])*(Σ[k,l] + μ[k]*μ[l]) - m([l])*μ̂[i]*μ[k] - m([k])*μ̂[i]*μ[l] - m0*μ̂[i]*(Σ[k,l] + μ[k]*μ[l]) for k in 1:n, l in 1:n]
    
    #before applying R_i
    Ib3(i) = [m([k,l]) - m([k])*μA[l] - m([l])*μA[k] + m0*(μA[k]*μA[l]-Σ̂[k,l]) for k in 1:n, l in 1:n]
    I3(i) = R(i)(Ib3(i))#QQQQ - do here - apply R_i


    function I4(i)
        mat = zeros(n,n)
        for l in 1:n
            for k in l:n #only k ≥ l
                ĩ = #QQQQ - fill in ĩĩ(k,l)
                j̃ = #QQQQ - fill in j̃j̃(k,l)
                #order 4 and 3
                mat[k,l] += m([j̃,ĩ,k,l]) - m([j̃,ĩ,k])*μ[l] - m([j̃,ĩ,l])*μ[k] - m([ĩ,k,l])*μA[j̃] - m([j̃,k,l])*μA[ĩ]

                #order 2
                mat[k,l] += m([j̃,ĩ])*(Σ[k,l] - μ[k]*μ[l]) + m([k,l])*(μA[j̃]*μA[ĩ] - Σ̂[j̃,ĩ]) - m([j̃,k])*μA[ĩ]*μ[l] - m([ĩ,l])*μA[j̃]*μ[k]

                #order 1 
                mat[k,l] += m([j̃])*μA[ĩ]*(μ[k]*μ[l] - Σ[k,l]) + m([ĩ])*μA[j̃]*(μ[k]*μ[l] - Σ[k,l]) + m([k])*μ[l]*(μA[j̃]*μA[ĩ] - Σ̂[j̃,ĩ]) + m([l])*μ[k]*(μA[j̃]*μA[ĩ] - Σ̂[j̃,ĩ])

                #order 0 
                mat[k,l] += m0*(μA[j̃]*μA[ĩ]*μ[k]*μ[l] - Σ̂[j̃,ĩ]*μ[k]*μ[l] - Σ[k,l]*μA[j̃]*μA[ĩ] + Σ̂[j̃, ĩ]*Σ[k,l])
            end
        end
        return mat
    end   

    return sum([I1(i) .* (I2(i)*U') for i in 1:n]) + sum([I3(i) .* (I4(i)*U') for i in 0:(n^2-1)])
end