const max_moment_levels = 2 #Just for mean and covariance matrix
const Children = Union{Vector{Vector{Int}},Nothing}

mutable struct BoxTruncatedMvNormal
    μₑ::Vector{Float64}     #mean parameter (of non-truncated form)
    Σₑ::Matrix{Float64}     #covariance parameter of non-truncated form)
    a::Vector{Float64}      #lower bounding box coordinates
    b::Vector{Float64}      #upper bounding box coordinates
    n::Int                #dimension

    #Children of type 'a' or 'b' are lower dimensional distributions used to for recursive computation
    children_a::Vector{BoxTruncatedMvNormal}
    children_b::Vector{BoxTruncatedMvNormal}

    #for each moment vector e.g. [0,1,0,1] or [0,0,2,0] has the tuple which is the computed (non-normalized) moment integral
    #of that vector and a list of children vectors
    momentDict::Dict{Vector{Int},Float64} #note that the values are non-normalized moment integrals
    treeDict::Dict{Vector{Int},Children}
    momentsComputed::Bool
    function BoxTruncatedMvNormal(  μₑ::Vector{Float64},
                                    Σₑ::Matrix{Float64}, 
                                    a::Vector{Float64}, 
                                    b::Vector{Float64}) 
        n = length(μₑ)
        length(a) != n && error("The length of a does not match the length of μₑ")
        length(b) != n && error("The length of b does not match the length of μₑ")
        size(Σₑ) != (n,n) && error("Dimension of Σₑ mismatch")
        !isposdef(Σₑ) && error("Σₑ is not positive definite")
        a > b && error("The a vector must be less than the b vector")

        if n ≥ 2
            μᵃ = [μₑ[setdiff(1:n,j)]  +  Σₑ[setdiff(1:n,j),j] * (a[j]-μₑ[j])/Σₑ[j,j] 
                    for j in 1:n]
            μᵇ = [μₑ[setdiff(1:n,j)]  +  Σₑ[setdiff(1:n,j),j] * (b[j]-μₑ[j])/Σₑ[j,j] 
                    for j in 1:n]
            Σ̃ = [Σₑ[setdiff(1:n,j),setdiff(1:n,j)] - (1/Σₑ[j,j])*Σₑ[setdiff(1:n,j),j]*Σₑ[j,setdiff(1:n,j)]' for j in 1:n]
            children_a = [BoxTruncatedMvNormal(μᵃ[j],Σ̃[j],a[setdiff(1:n,j)],b[setdiff(1:n,j)]) for j in 1:n]
            children_b = [BoxTruncatedMvNormal(μᵇ[j],Σ̃[j],a[setdiff(1:n,j)],b[setdiff(1:n,j)]) for j in 1:n]
        else #n==1
            children_a = Array{BoxTruncatedMvNormal,1}[] #no children
            children_b = Array{BoxTruncatedMvNormal,1}[] #no children
        end
        momentDict, treeDict = init_dicts(n)
        new(μₑ,Σₑ,a,b,n,children_a,children_b,momentDict,treeDict,false)
    end
end

function Base.show(io::IO, d::BoxTruncatedMvNormal) 
    println(io, "Box Truncated MvNormal")
    println(io, "n = $(d.n)")
    println(io, "μₑ = $(d.μₑ)" )
    println(io, "Σₑ = $(d.Σₑ)" )
    println(io, "α = $(alpha(d))")
    println("Limits:")
    
    for i in 1:d.n
        println(io, "$i:\t ",(d.a[i],d.b[i]))
    end
    
    if d.momentsComputed
        println("Moments:")
        for k in keys(d.momentDict)
            println(io, k, "\t", moment(d,k))
        end
        println(io,"mean:", mean(d))
        println(io,"cov:", cov(d))
    else
        println("Moments not computed")
    end
end

function init_dicts(n)
    function addToBaseKey(  baseKey::Vector{Int},
                            n::Int,
                            md::Dict{Vector{Int},Float64},
                            td::Dict{Vector{Int},Children})
        keys = Vector{Vector{Int}}(undef,n)
        for i in 1:n
            key = copy(baseKey)
            key[i] += 1
            md[key] = NaN
            td[key] = nothing
            keys[i] = key
        end
        md[baseKey] = NaN
        td[baseKey] = keys
        keys
    end

    md = Dict{Vector{Int},Float64}() #momentDict
    td = Dict{Vector{Int},Children}() #treeDict
    rootKey = zeros(Int,n)
    key_vals = [rootKey]
    for _ = 1:max_moment_levels
        levelKeys = Vector{Int}[]
        for key in key_vals
            newKeys = addToBaseKey(key,n,md,td)
            append!(levelKeys,newKeys)
        end
        key_vals = levelKeys
    end
    md,td
end

function compute_moments(d::BoxTruncatedMvNormal)
    function compute_children_moments(d::BoxTruncatedMvNormal,baseKey::Vector{Int})
        d.treeDict[baseKey] == nothing && return #recursion stopping criteria
        c = c_vector(d,baseKey)
        for k in d.treeDict[baseKey]
            i = findfirst((x)->x==1,k-baseKey)     
            d.momentDict[k] = d.μₑ[i]*d.momentDict[baseKey] + (d.Σₑ*c)[i]
            compute_children_moments(d,k)
        end
    end

    function c_vector(d::BoxTruncatedMvNormal,k::Vector{Int})
        c = Vector{Float64}(undef,d.n)
        for j in 1:d.n
            kMinus = copy(k)
            kMinus[j] = kMinus[j]-1
            F0 = (sum(kMinus .>= 0) == d.n) ? d.momentDict[kMinus] : 0.0 
            F1 = length(d.children_a) >= 1 ? raw_moment(d.children_a[j],k[setdiff(1:d.n,j)]) : 0.0
            F2 = length(d.children_b) >= 1 ? raw_moment(d.children_b[j],k[setdiff(1:d.n,j)]) : 0.0
            ϕ = pdf.(Normal(d.μₑ[j],sqrt(d.Σₑ[j,j])),[d.a[j],d.b[j]] )
            c[j] =  k[j]*F0  +  d.a[j]^k[j]*ϕ[1]*F1  -  d.b[j]^k[j]*ϕ[2]*F2 
        end
        c
    end
    
    if d.n > 1
        baseKey = zeros(Int,d.n) #[0,0,....,0]
        d.momentDict[baseKey] = LL(d)
        compute_children_moments(d,baseKey) #start recursion
    else  #n==1
        #TODO - Improve this code to go beyond max_moment_levels =2
        @assert max_moment_levels == 2
        @assert d.n == 1
        distTruncated = TruncatedNormal(d.μₑ[1],sqrt(d.Σₑ[1]),d.a[1],d.b[1])
        distFull = Normal(d.μₑ[1],sqrt(d.Σₑ[1]))
        #QQQQ Replace with recursive one dimensional
        d.momentDict[[0]] = cdf(distFull,d.b[1])-cdf(distFull,d.a[1])
        d.momentDict[[1]] = Distributions.mean(distTruncated)*d.momentDict[[0]]
        d.momentDict[[2]] = (var(distTruncated) + Distributions.mean(distTruncated)^2)*d.momentDict[[0]]
    end
    d.momentsComputed = true
end

function raw_moment(d::BoxTruncatedMvNormal,k::Vector{Int})
    !d.momentsComputed && compute_moments(d)
    return d.momentDict[k]
end

moment(d::BoxTruncatedMvNormal,k::Vector{Int}) = raw_moment(d,k) / alpha(d)

alpha(d::BoxTruncatedMvNormal) = raw_moment(d,zeros(Int,d.n))

function mean(d::BoxTruncatedMvNormal)
    μ = Vector{Float64}(undef,d.n)
    for i in 1:d.n
        ee = zeros(Int,d.n)
        ee[i] = 1
        μ[i] = moment(d,ee)
    end
    μ
end

function cov(d::BoxTruncatedMvNormal)
    Σ = zeros(Float64,d.n,d.n)
    for i in 1:d.n, j in 1:d.n
        ee = zeros(Int,d.n)
        if i == j
            ee[i] = 2
        else
            ee[i], ee[j] = 1, 1
        end
        Σ[i,j] = moment(d,ee)
    end
    μ = mean(d)
    Σ-μ*μ'
end

function rand(d::BoxTruncatedMvNormal)
    rand(MvNormal(d.μₑ,d.Σₑ)) #TODO QQQQ
end


function pdf_nontruncated(d::BoxTruncatedMvNormal,x)
    d_nontruncated = MvNormal(d.μₑ,d.Σₑ)
    pdf(d_nontruncated,x)
end