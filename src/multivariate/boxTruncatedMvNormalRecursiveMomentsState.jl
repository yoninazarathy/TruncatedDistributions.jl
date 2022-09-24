# const max_moment_levels = 2 #Just for mean and covariance matrix
const Children = Union{Vector{Vector{Int}},Nothing}

mutable struct BoxTruncatedMvNormalRecursiveMomentsState <: TruncatedMvDistributionState
    d::MvNormal
    r::BoxTruncationRegion
    n::Int                  #dimension
    max_moment_levels::Int
    
    # # #Children of type 'a' or 'b' are lower dimensional distributions used to for recursive computation
    children_a::Vector{BoxTruncatedMvNormalRecursiveMomentsState}
    children_b::Vector{BoxTruncatedMvNormalRecursiveMomentsState}

    # # #for each moment vector e.g. [0,1,0,1] or [0,0,2,0] has the tuple which is the computed (non-normalized) moment integral
    # # #of that vector and a list of children vectors
    rawMomentDict::Dict{Vector{Int},Float64} #note that the values are non-normalized moment integrals
    treeDict::Dict{Vector{Int},Children}
    rawMomentsComputed::Bool

    tp::Float64
    μ::Vector{Float64}
    Σ::PDMat
    tp_err::Float64
    μ_err::Float64
    Σ_err::Float64

    function BoxTruncatedMvNormalRecursiveMomentsState(d::MvNormal, r::BoxTruncationRegion, max_moment_levels::Int)
        μₑ, Σₑ = d.μ, d.Σ
        a, b = r.a, r.b                              
        n = length(d)
        length(a) != n && error("The length of a does not match the length")
        length(b) != n && error("The length of b does not match the length")
        a > b && error("The a vector must be less than the b vector")

        if n ≥ 2
            μᵃ = [μₑ[setdiff(1:n,j)]  +  Σₑ[setdiff(1:n,j),j] * (a[j]-μₑ[j])/Σₑ[j,j] 
                    for j in 1:n]
            μᵇ = [μₑ[setdiff(1:n,j)]  +  Σₑ[setdiff(1:n,j),j] * (b[j]-μₑ[j])/Σₑ[j,j] 
                    for j in 1:n]
            Σ̃ = [Σₑ[setdiff(1:n,j),setdiff(1:n,j)] - (1/Σₑ[j,j])*Σₑ[setdiff(1:n,j),j]*Σₑ[j,setdiff(1:n,j)]' for j in 1:n]
            children_a = [BoxTruncatedMvNormalRecursiveMomentsState(
                            MvNormal(μᵃ[j],Σ̃[j]),
                            BoxTruncationRegion(a[setdiff(1:n,j)],b[setdiff(1:n,j)]),
                            max_moment_levels) for j in 1:n]
            children_b = [BoxTruncatedMvNormalRecursiveMomentsState(
                            MvNormal(μᵇ[j],Σ̃[j]),
                            BoxTruncationRegion(a[setdiff(1:n,j)],b[setdiff(1:n,j)]),
                            max_moment_levels) for j in 1:n]
        else #n==1
            children_a = Array{BoxTruncatedMvNormalRecursiveMomentsState,1}[] #no children
            children_b = Array{BoxTruncatedMvNormalRecursiveMomentsState,1}[] #no children
        end
        rawMomentDict, treeDict = init_dicts(n,max_moment_levels)
        new(d,r,n,max_moment_levels,children_a,children_b,rawMomentDict,treeDict,false,
            NaN,
            Vector{Float64}(undef,0),
            PDMat(Array{Float64,2}(I,n,n)),
            Inf, Inf, Inf)
    end
end

function init_dicts(n::Int,max_moment_levels::Int)
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

    md = Dict{Vector{Int},Float64}() #rawMomentDict
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

function compute_moments(d::BoxTruncatedMvNormalRecursiveMomentsState)
    function compute_children_moments(d::BoxTruncatedMvNormalRecursiveMomentsState,baseKey::Vector{Int})
        d.treeDict[baseKey] == nothing && return #recursion stopping criteria
        c = c_vector(d,baseKey)
        for k in d.treeDict[baseKey]
            i = findfirst((x)->x==1,k-baseKey)     
            d.rawMomentDict[k] = d.d.μ[i]*d.rawMomentDict[baseKey] + (d.d.Σ*c)[i]
            compute_children_moments(d,k)
        end
    end

    function c_vector(d::BoxTruncatedMvNormalRecursiveMomentsState,k::Vector{Int})
        c = Vector{Float64}(undef,d.n)
        for j in 1:d.n
            kMinus = copy(k)
            kMinus[j] = kMinus[j]-1
            F0 = (sum(kMinus .>= 0) == d.n) ? d.rawMomentDict[kMinus] : 0.0 
            F1 = length(d.children_a) >= 1 ? raw_moment(d.children_a[j],k[setdiff(1:d.n,j)]) : 0.0
            F2 = length(d.children_b) >= 1 ? raw_moment(d.children_b[j],k[setdiff(1:d.n,j)]) : 0.0
            ϕ = pdf.(Normal(d.d.μ[j],sqrt(d.d.Σ[j,j])),[d.r.a[j],d.r.b[j]] )
            c[j] =  k[j]*F0  +  d.r.a[j]^k[j]*ϕ[1]*F1  -  d.r.b[j]^k[j]*ϕ[2]*F2 
        end
        c
    end
    
    if d.n > 1
        baseKey = zeros(Int,d.n) #[0,0,....,0]
        d.rawMomentDict[baseKey] = LL(d)
        compute_children_moments(d,baseKey) #start recursion
    else  #n==1
        @assert d.n == 1
        distTruncated = TruncatedNormal(d.d.μ[1],sqrt(d.d.Σ[1]),d.r.a[1],d.r.b[1])
        d.rawMomentDict[[0]] = distTruncated.tp
        m = moments(distTruncated, d.max_moment_levels)
        for i in 1:d.max_moment_levels
            d.rawMomentDict[[i]] = m[i]*distTruncated.tp
        end
    end
    d.rawMomentsComputed = true
end

function raw_moment(d::BoxTruncatedMvNormalRecursiveMomentsState,k::Vector{Int})
    !d.rawMomentsComputed && compute_moments(d)
    return d.rawMomentDict[k]
end

function raw_moment_dict(d::BoxTruncatedMvNormalRecursiveMomentsState)
    !d.rawMomentsComputed && compute_moments(d)
    return copy(d.rawMomentDict)
end


# moment(d::BoxTruncatedMvNormalRecursiveMomentsState,k::Vector{Int}) = raw_moment(d,k) / alpha(d)

# alpha(d::BoxTruncatedMvNormalRecursiveMomentsState) = raw_moment(d,zeros(Int,d.n))

# function mean(d::BoxTruncatedMvNormalRecursiveMomentsState)
#     μ = Vector{Float64}(undef,d.n)
#     for i in 1:d.n
#         ee = zeros(Int,d.n)
#         ee[i] = 1
#         μ[i] = moment(d,ee)
#     end
#     μ
# end

# function cov(d::BoxTruncatedMvNormalRecursiveMomentsState)
#     Σ = zeros(Float64,d.n,d.n)
#     for i in 1:d.n, j in 1:d.n
#         ee = zeros(Int,d.n)
#         if i == j
#             ee[i] = 2
#         else
#             ee[i], ee[j] = 1, 1
#         end
#         Σ[i,j] = moment(d,ee)
#     end
#     μ = mean(d)
#     Σ-μ*μ'
# end

# function rand(d::BoxTruncatedMvNormalRecursiveMomentsState)
#     rand(MvNormal(d.μₑ,d.Σₑ)) #TODO QQQQ
# end


# function pdf_nontruncated(d::BoxTruncatedMvNormalRecursiveMomentsState,x)
#     d_nontruncated = MvNormal(d.μₑ,d.Σₑ)
#     pdf(d_nontruncated,x)
# end

function LL(d::BoxTruncatedMvNormalRecursiveMomentsState)
    @info "doing base numerical integral on dimension $(d.n)."
    hcubature((x)->pdf(d.d,x),d.r.a,d.r.b,maxevals = 10^6)[1]
end