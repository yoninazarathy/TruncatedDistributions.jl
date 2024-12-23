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
    I3 = [(
        m([ii(k), jj(k)]) - m([ii(k)])*μA[jj(k)] - m([jj(k)])*μA[ii(k)] 
        + m0*(μA[ii(k)]*μA[jj(k)] - Σ̂[ii(k),jj(k)])
        )
        for k in 1:n^2]
    I4 = [(
        m([ii(k),jj(k),l]) - m([ii(k),l])*μA[jj(k)] - m([jj(k),l])*μA[ii(k)] - m([ii(k),jj(k)])*μ[l] 
        + m([l])*(μA[ii(k)]*μA[jj(k)]-Σ̂[ii(k),jj(k)]) + m([jj(k)])*μA[ii(k)]*μ[l] 
        - m0*(μA[ii(k)]*μA[jj(k)] - Σ̂[ii(k),jj(k)]) 
        ) 
        for k in 1:n^2, l in 1:n]
    return ((I1'*I2 + I3'*I4)*inv(Matrix(Σ)))'
end

function U_gradient(    d::RecursiveMomentsBoxTruncatedMvNormal, 
                        μA::Vector{Float64}, 
                        μ̂::Vector{Float64}, 
                        Σ̂::PDMat;
                        U::Union{Nothing, Matrix{Float64}} = nothing) # U is optional
    μ = d.untruncated.μ
    Σ = d.untruncated.Σ
    if isnothing(U)
        U = cholesky(0.5*(inv(Matrix(Σ)) + inv(Matrix(Σ))')).U
    end
    n = length(d)

    m(inds) = raw_moment_from_indices(d, inds)
    m0 = m(Int[])
    # testm0 = hcubature((x)->pdf(d.untruncated,x),d.region.a, d.region.b; atol=0.00001)[1]
    # @show m0-testm0

    # m1 = m([1])
    # testm1 = hcubature((x)->pdf(d.untruncated,x)*x[1],d.region.a, d.region.b; atol=0.00001)[1]
    # @show m1-testm1

    # m2 = m([2])
    # testm2 = hcubature((x)->pdf(d.untruncated,x)*x[2],d.region.a, d.region.b; atol=0.00001)[1]
    # @show m2-testm2

    # m12 = m([1,2])
    # testm12 = hcubature((x)->pdf(d.untruncated,x)*x[1]*x[2],d.region.a, d.region.b; atol=0.00001)[1]
    # @show m12-testm12

    # m11 = m([1,1])
    # testm11 = hcubature((x)->pdf(d.untruncated,x)*x[1]*x[1],d.region.a, d.region.b; atol=0.00001)[1]
    # @show m11-testm11

    # m22 = m([2,2])
    # testm22 = hcubature((x)->pdf(d.untruncated,x)*x[2]*x[2],d.region.a, d.region.b; atol=0.00001)[1]
    # @show m22-testm22

    # m112 = m([1,1,2])
    # testm112 = hcubature((x)->pdf(d.untruncated,x)*x[1]*x[1]*x[2],d.region.a, d.region.b; atol=0.00001)[1]
    # @show m112-testm112

    # m1122 = m([1,1,2,2])
    # testm1122 = hcubature((x)->pdf(d.untruncated,x)*x[1]*x[1]*x[2]*x[2],d.region.a, d.region.b; atol=0.00001)[1]
    # @show m1122-testm1122

    I1(i) = fill(m([i]) - m0*μ̂[i], n, n)
    # CODE for testing I1
    # I1test(i) = hcubature((x)->pdf(d.untruncated,x)*fill(x[i] - μ̂[i], n, n),d.region.a, d.region.b; atol=0.00001)[1]
    # @show maximum(abs, [norm(I1(1)- I1test(1)) for i in 1:n])
    # I2(i) = [k < l ? 0 : (
    #             -m([i,k,l]) 
    #             + m([i,k])*μ[l] + m([k,l])*μ̂[i] + m([i,l])*μ[k] 
    #             - m([i])*(Σ[k,l] + μ[k]*μ[l]) - m([l])*μ̂[i]*μ[k] - m([k])*μ̂[i]*μ[l] 
    #             + m0*μ̂[i]*(Σ[k,l] + μ[k]*μ[l])
    #             ) for k in 1:n, l in 1:n]
    # CODE for testing I2
    Q = triu(fill(1,n,n))
    iU = inv(U)
    # @show iU
    # @show Σ
    # @show iU*iU'
    # I2test(i) = (hcubature((x)->pdf(d.untruncated,x)*fill(μ̂[i] - x[i], n, n).*(Q .* (iU' + (U*(x-μ))*(x-μ)')),d.region.a, d.region.b; atol=0.0000001)[1])
    # I2trans(i) = (hcubature((x)->pdf(d.untruncated,x)*fill(μ̂[i] - x[i], n, n).*(Q .* (iU' + (U*(x-μ))*(x-μ)')'),d.region.a, d.region.b; atol=0.0000001)[1])
    # println("I2(1):");display(I2(1))
    # println("I2trans(1):");display(I2trans(1))
    # I2b(i) = (hcubature((x)->pdf(d.untruncated,x)*fill(μ̂[i] - x[i], n, n) .* Q' .* (iU + (x-μ)*(x-μ)'*U'),d.region.a, d.region.b; atol=0.0000001)[1])
    # I2c(i) = (hcubature((x)->pdf(d.untruncated,x)*fill(μ̂[i] - x[i], n, n) .* Diagonal(iU),d.region.a, d.region.b; atol=0.0000001)[1]) - Q' .* (hcubature((x)->pdf(d.untruncated,x)*fill(x[i]- μ̂[i], n, n) .* ((x-μ)*(x-μ)'*U'),d.region.a, d.region.b; atol=0.0000001)[1])
    # @show iU

    I2(i) = [k > l ? 0.0 : (k == l ? iU[k,k]*(m0*μ̂[i] - m([i])) : 0.0) - sum(U[k,j]*(
        m([i,l,j]) - m([i,j])*μ[l] - m([l,j])*μ̂[i] - m([i,l])*μ[j] + m([i])*μ[l]*μ[j] + m([j])*μ[l]*μ̂[i] + m([l])*μ̂[i]*μ[j] - m0*μ̂[i]*μ[l]*μ[j]
        )  for j in k:n) for k=1:n, l=1:n]
    
    
    # println("I2(1):");display(I2(1))
    # println("Itest(1):");display(I2test(1))
    # println("I2c(2):");display(I2c(2))
    # println("I2new(2):");display(I2new(2))
    # I2testb(i) = (hcubature((x)->pdf(d.untruncated,x)*fill(μ̂[i] - x[i], n, n).* Q'.* (Σ + (x-μ)*(x-μ)'),d.region.a, d.region.b; atol=0.0000001)[1])*U'
    # println("I2(1):");display(I2(1))
    # println("I2(1):");display(I2(1))
    # println("I2test(1):");display(I2test(1))
    # @show maximum(abs, [norm(I2(i)- I2test(i)) for i in 1:n])

    #before applying R_i
    Ib3 = [m([k,l]) - m([k])*μA[l] - m([l])*μA[k] + m0*(μA[k]*μA[l]-Σ̂[k,l]) for k in 1:n, l in 1:n]
    
    Rby1(A) = reshape(vcat(vec(A)[end],vec(A)[1:end-1]),size(A)...)
    function R(A,i)
        B = copy(A)
        for _ in 1:i
            B = Rby1(B)
        end
        return B
    end
    # I3(i) = R(Ib3,i)

    shiftindex(k,l,i,n) = mod((l-1)*n+(k-1)-i, n^2)
    itild(k,l,i,n) = mod(shiftindex(k,l,i,n),n)+1
    jtild(k,l,i,n) = div(shiftindex(k,l,i,n),n)+1

    # A = [10 30;
    #      20 40]

    # # A = [10 40 70;
    # #      20 50 80;
    # #     30  60 90]
    # print("R_0:");display(R(A,0))
    # print("R_1:");display(R(A,1))
    # print("R_2:");display(R(A,2))
    # print("R_3:");display(R(A,3))
    # print("R_4:");display(R(A,4))

    # for i in 0:(n^2-1)
    #     @show i
    #     for l in 1:n
    #         for k in 1:n
    #             si = shiftindex(k,l,i,n)
    #             ii, jj = itild(k,l,i,n), jtild(k,l,i,n)
    #             @show (k,l) => si,(ii,jj) 
    #         end
    #     end
    # end

    # function I4(i)
    #     mat = zeros(n,n)
    #     for l in 1:n
    #         for k in l:n #only k ≥ l
    #             ĩ = itild(k,l,i,n)
    #             j̃ = jtild(k,l,i,n)

    #             #order 4 and 3
    #             mat[k,l] += -m([j̃,ĩ,k,l]) + m([j̃,ĩ,k])*μ[l] + m([j̃,ĩ,l])*μ[k] + m([ĩ,k,l])*μA[j̃] + m([j̃,k,l])*μA[ĩ]

    #             #order 2
    #             mat[k,l] += -m([j̃,ĩ])*(Σ[k,l] + μ[k]*μ[l]) + m([k,l])*(Σ̂[j̃,ĩ]- μA[j̃]*μA[ĩ]) - m([j̃,k])*μA[ĩ]*μ[l] - m([ĩ,l])*μA[j̃]*μ[k] - m([j̃,l])*μA[ĩ]*μ[k] - m([ĩ,k])*μA[j̃]*μ[l]

    #             #order 1 
    #             mat[k,l] += m([j̃])*μA[ĩ]*(Σ[k,l]+μ[k]*μ[l]) + m([ĩ])*μA[j̃]*(Σ[k,l] + μ[k]*μ[l]) - m([k])*μ[l]*(Σ̂[j̃,ĩ] - μA[j̃]*μA[ĩ]) - m([l])*μ[k]*(Σ̂[j̃,ĩ] - μA[j̃]*μA[ĩ])

    #             #order 0 
    #             mat[k,l] += m0*(Σ̂[j̃,ĩ] - μA[j̃]*μA[ĩ])*(Σ[k,l] -μ[k]*μ[l])
    #         end
    #     end
    #     return mat
    # end   

    # I4(i) = (hcubature((x)->pdf(d.untruncated,x)*R((x-μA)*(x-μA)'-Σ̂,i).*((-Q) .* (iU' + (U*(x-μ))*(x-μ)')),d.region.a, d.region.b; atol=0.000000001)[1])

    # I4num(i) = hcubature((x)->pdf(d.untruncated,x) * Q .* R((x-μA)*(x-μA)'-Σ̂,i) .* (U*(x-μ)*(x-μ)'),d.region.a, d.region.b; maxevals = 10^6)[1]#atol=0.00001)[1]

    function MM1(i,k,j,l)
        ans = 0.0
        ĩ = itild(k,l,i,n)
        j̃ = jtild(k,l,i,n)

        #order 4
        ans += m([j̃,ĩ,l,j]) 

        #order 3
        ans -=  m([j̃,ĩ,l])*μ[j] + m([j̃,ĩ,j])*μ[l] + m([ĩ,l,j])*μA[j̃] + m([j̃,l,j])*μA[ĩ]

        #order 2
        ans += m([j̃,ĩ])*μ[l]*μ[j] + m([l,j̃])*μA[ĩ]*μ[j] + m([ĩ,j])*μA[j̃]*μ[l] + m([j̃,j])*μA[ĩ]*μ[l] + m([ĩ,l])*μA[j̃]*μ[j] + m([j,l])*μA[j̃]*μA[ĩ]

        #order 1 
        ans -= m([j̃])*μA[ĩ]*μ[l]*μ[j] + m([ĩ])*μA[j̃]*μ[l]*μ[j] + m([l])*μ[j]*μA[j̃]*μA[ĩ] + m([j])*μ[l]*μA[j̃]*μA[ĩ]

        #order 0 
        ans += m0*μA[j̃]*μA[ĩ]*μ[l]*μ[j]

        return ans
    end   

    function MM2(i,k,j,l)
        ĩ = itild(k,l,i,n)
        j̃ = jtild(k,l,i,n)

        return Σ̂[ĩ,j̃]*(m([l,j]) - m([l])*μ[j] - m([j])*μ[l] + m0*μ[l]*μ[j])
    end

    I4(i) = [k <= l ? sum(U[k,j]*(MM1(i,k,j,l) - MM2(i,k,j,l)) for j in 1:n) : 0  for k in 1:n, l in 1:n]
# 
    #  I4(i) = hcubature((x)->pdf(d.untruncated,x) * Q .* R((x-μA)*(x-μA)'-Σ̂,i) .* (U*(x-μ)*(x-μ)'),d.region.a, d.region.b; maxevals = 10^6)[1]#atol=0.00001)[1]


    # I4(i) = (-Diagonal(iU) .* I3(i)' - missing4_test(i))'
    # println("I4(2):");display(I4(2))
    # println("I4num(2):");display(I4num(2))

    return sum([I1(i) .* I2(i) for i in 1:n]) - sum([R(Ib3,i) .* (Diagonal(iU) .* R(Ib3,i) + I4(i)) for i in 0:(n^2-1)])
    # return triu(sum([I1(i) .* I2(i)*(U')) for i in 1:n]) + sum([I3(i) .* (U*(I4(i)')) for i in 0:(n^2-1)]))
end