# TruncatedDistributions.jl Documentation

A Julia package for Truncated Distributions including multi-variate truncated distributions.

For example:
```jldoctest
julia> using TruncatedDistributions, PDMats

julia> d = BasicBoxTruncatedMvNormal([2.3,4.3],PDMat([1.0 0.3; 0.3  2.4]),[-1.0, 0.5],[5.2,1.8]);

julia> mean(d)
2-element Vector{Float64}:
 1.930922186217924
 1.3227326677503652

julia> d.state.μ_err
3.457541323608117e-8
```


```@contents
```


## Specific Types

```@docs
BasicBoxTruncatedMvNormal
RecursiveMomentsBoxTruncatedMvNormal
```

## Common Types

The key type exposed by this package is `TruncatedMvDistribution`. It is parameterized by a `MultivariateDistribution` from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl), a `TruncationRegion`, and an additional type called the `TruncatedMvDistributionState`.

```@docs
TruncatedMvDistribution
TruncationRegion
TruncatedMvDistributionState
```

## Functions

```@docs
TruncatedMvDistribution{D <: MultivariateDistribution, R <: TruncationRegion}
```

## Index

```@index
```

