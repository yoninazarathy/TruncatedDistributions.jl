function normal2a()
    μₑ = [2.5, 3.5]
    Σₑ = [2.0 -0.5;
        -0.5 5.0]
    a = [-2.3,-20.]
    b = [4.,12.3]
    properties = Dict{String,Any}()
    properties["length"] = 2
    properties["tp"] = 0.8551938607791414
    properties["mean"] = [2.126229598541412, 3.5930224468784577]
    properties["covariance"] = [1.0 0; 0 1.0] #QQQQ
    properties["moment123"] = 2.3 #QQQQ
    (MvNormal(μₑ,Σₑ), BoxTruncationRegion(a,b),properties)
end

function normal2b()
    μₑ = [2.5, 3.5]
    Σₑ = [  3.3 0.5;
            0.5 5.0]
    a = [-5.4,-20.]
    b = [2.4,6.3]
    properties = Dict{String,Any}()
    properties["length"] = 2
    properties["tp"] = 0.43660920327458974
    properties["mean"] = [0.9734400003512856, 2.886032492774952]
    properties["covariance"] = [1.0 0; 0 1.0] #QQQQ
    properties["moment123"] = 2.3 #QQQQ
    (MvNormal(μₑ,Σₑ), BoxTruncationRegion(a,b),properties)
end

function normal2c()
    a = [-1.5, -1.4]
    b = [1.7, 2.4]
    μ₁ = 0.45; μ₂ = -0.2
    σ₁ = 1.2; σ₂ = 0.8
    ρ = 0.6
    μₑ = [μ₁, μ₂]
    Σₑ = [σ₁^2 ρ*σ₁*σ₂ ;ρ*σ₁*σ₂ σ₂^2];
    properties = Dict{String,Any}()
    properties["length"] = 2
    properties["tp"] = 0.7517763397386328
    properties["mean"] = [0.308275927584875, -0.18541888959459515]
    properties["covariance"] = [1.0 0; 0 1.0] #QQQQ
    properties["moment123"] = 2.3 #QQQQ
    (MvNormal(μₑ,Σₑ), BoxTruncationRegion(a,b),properties)
end

function normal3a()
    μₑ = [3.5,2,3.5]
    Σₑ = [ 7. 1 0 ; 
           1 3.3 2  ;
           0 2 3.8 ]
    a = [-2. ,-1 ,2]
    b = [4. ,4. , 5]
    properties = Dict{String,Any}()
    properties["length"] = 3
    properties["tp"] = 0.9 #QQQQ
    properties["mean"] = [1.0,1.0] #QQQQ
    properties["covariance"] = [1.0 0; 0 1.0] #QQQQ
    properties["moment123"] = 2.3 #QQQQ
    (MvNormal(μₑ,Σₑ), BoxTruncationRegion(a,b),properties)
end

function normal4a()
    μₑ = [3.5,2,3.5,3.5]
    Σₑ = [ 7. 1 0 1 ; 
        1 3.3 2 0 ;
        0 2 3.4 0 ;
        1 0 0 4 ]
    a = [-20. ,-20, -20 ,-20]
    b = [20. ,20. ,20, 20]
    properties = Dict{String,Any}()
    properties["length"] = 4
    properties["tp"] = 0.9 #QQQQ
    properties["mean"] = [1.0,1.0] #QQQQ
    properties["covariance"] = [1.0 0; 0 1.0] #QQQQ
    properties["moment123"] = 2.3 #QQQQ
    (MvNormal(μₑ,Σₑ), BoxTruncationRegion(a,b),properties)
end

function normal5a()
    μₑ = [3.5,2,3.5,3.5,5.3]
    Σₑ = [ 7. 1   0 0 1 ; 
        1  3.3 2 0 0 ;
        0 2 3.4 0 0 ;
        0 0 0 4 0; 
        1 0 0 0 2]
    a = [-20.,-20,-20,-20,-20]
    b = [20.,20.,20,20,20]
    properties = Dict{String,Any}()
    properties["length"] = 5
    properties["tp"] = 0.9 #QQQQ
    properties["mean"] = [1.0,1.0] #QQQQ
    properties["covariance"] = [1.0 0; 0 1.0] #QQQQ
    properties["moment123"] = 2.3 #QQQQ
    (MvNormal(μₑ,Σₑ), BoxTruncationRegion(a,b),properties)
end

distribution_generators2 = [normal2a,normal2b,normal2c]
distribution_generators3 = [normal3a]
distribution_generators4 = [normal4a]
distribution_generators5 = [normal5a]
distribution_generators = vcat(distribution_generators2,distribution_generators3,distribution_generators4,distribution_generators5)