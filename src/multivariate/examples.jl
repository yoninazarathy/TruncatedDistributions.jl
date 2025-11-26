"""
    Returns an example.
"""
function get_example(;dist_type = MvNormal, n = 2, index = 1)
    if dist_type == MvNormal
        if haskey(normal_examples, n)
            return normal_examples[n][index]
        else
            error("No examples for $dist_type with n=$n")
        end
    else
        error("Not supporting dist type $dist_type")
    end
end

"""
    Returns the number of examples for a given n.
"""
function get_num_examples(n; dist_type = MvNormal)
    if dist_type == MvNormal
        return length(normal_examples[n])
    else
        error("Not supporting dist type $dist_type")
    end
end

"""
    Returns the possible example sizes.
"""
function get_example_sizes(;dist_type = MvNormal)
    return sort(collect(keys(normal_examples)))
end

"""
    An example. The `arb_moment_to_check_index` is some tuple of indexes of a moment to check (can be nothing). 
    The associated value is `arb_moment_to_check_value`. This is for checking some single (arbitrary) moment.
"""
@with_kw struct NormalExample
    n::Int = 2
    μ::AbstractVector{Float64}
    Σ::AbstractMatrix{Float64}
    a::AbstractVector{Float64}
    b::AbstractVector{Float64}
    tp::Float64
    μ̂::AbstractVector{Float64}
    Σ̂::AbstractMatrix{Float64}
    arb_moment_to_check_index::Union{Nothing,Tuple{Vararg{Int}}} = nothing
    arb_moment_to_check_value::Union{Nothing, Float64} = nothing
end

"""
Create a distribution from an example.
"""
function dist_from_example(ne::NormalExample)
    return RecursiveMomentsBoxTruncatedMvNormal(ne.μ, PDMat(ne.Σ), ne.a, ne.b)
end

#########################
## Hard coded examples ##
#########################
normal_examples = Dict()
normal_examples[2] = [
    NormalExample(  μ = [2.5, 3.5],
                    Σ = [2.0 -0.5;
                        -0.5 5.0],
                    a = [-2.3,-20.],
                    b = [4.,12.3],
                    tp = 0.8551938607791414,
                    μ̂ = [2.126229598541412, 3.5930224468784577],
                    Σ̂ = [1.0 0; 0 1.0],
                    arb_moment_to_check_index = (3, 5),
                    arb_moment_to_check_value = 58529.1327440061),
    NormalExample(  μ = [2.5, 3.5],
                    Σ = [  3.3 0.5;
                            0.5 5.0],
                    a = [-5.4,-20.],
                    b = [2.4,6.3],
                    tp = 0.43660920327458974,
                    μ̂ = [0.9734400003512856, 2.886032492774952],
                    Σ̂ = [1.0 0; 0 1.0],
                    arb_moment_to_check_index = (0, 3),
                    arb_moment_to_check_value = 52.44849808917711)
]
normal_examples[3] = [
    NormalExample(  μ = [3.5,2,3.5],
                    Σ = [ 7. 1 0 ; 
                          1 3.3 2  ;
                          0 2 3.8 ],
                    a = [-4. ,-3 ,-1],
                    b = [7.5 ,6.5 , 6.5],
                    tp = 0.2771862142891479,
                    μ̂ = [3.1593375223480122, 1.8453845525318782, 3.3023816830081723],
                    Σ̂ =  [5.28031    0.719954  -0.0155715 ; 0.719954   2.8325     1.38021 ;-0.0155715  1.38021    2.71015],
                    arb_moment_to_check_index = (1, 1, 1),
                    arb_moment_to_check_value = 11.740894033054031)
]
normal_examples[4] = [
    NormalExample(  μ = [3.5,2,3.5,3.5],
                    Σ = [ 7. 1 0 1 ; 
                            1 3.3 2 0 ;
                            0 2 3.4 0 ;
                            1 0 0 4 ],
                    a = [-8.0 ,-10, -4 ,-20],
                    b = [6.0 ,15. ,13, 10],
                    tp = 0.2806670136910537,
                    μ̂ = [2.2496049084043497, 0.7606400499117575, 1.6967554685821928, 3.321372131126033],
                    Σ̂ = [1.0 0 0 0 ; 0 1.0 0 0; 0 0 1.0 0; 0 0 0 1.0],
                    arb_moment_to_check_index = (2, 2, 2, 3),
                    arb_moment_to_check_value = 10550.695322644422)
]
normal_examples[5] = []
normal_examples[10] = []
normal_examples[20] = []
normal_examples[50] = []
