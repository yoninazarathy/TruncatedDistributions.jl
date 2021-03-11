# TruncatedDistributions

This package provides support for univariate and multivariate distributions in Julia beyond the [univariate truncated distributions support](https://juliastats.org/Distributions.jl/latest/truncate/) that is in [Distributions.jl](https://github.com/JuliaStats/Distributions.jl), Julia's main distributions package. Beyond basic distribution functionality, the package provides functions for fitting parameters to desired moments.

At the moment the only supported multivariate distribution is Multivariate Normal with box truncation (truncation over a region). For univariate distributions, essentially any distribution can be truncated using the truncation mechanism from Distributions.jl and the functions of this package can be applied.  