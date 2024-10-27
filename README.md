# vfi.jl

<!-- Tidyverse lifecycle badges, see https://www.tidyverse.org/lifecycle/ Uncomment or delete as needed. -->
![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
<!-- 
![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) 
[![build](https://github.com/avinnofaruk/vfi.jl/workflows/CI/badge.svg)](https://github.com/avinnofaruk/vfi.jl/actions?query=workflow%3ACI)
-->
<!-- Documentation -- uncomment or delete as needed -->
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://avinnofaruk.github.io/vfi.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://avinnofaruk.github.io/vfi.jl/dev)
-->
<!-- Aqua badge, see test/runtests.jl -->
<!-- [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) -->

## About

Generic solver for the optimal neoclassical growth model, using the value function iteration (VFI) method. Currently accepts single factor of production $K$ (i.e. $L=1$), and single output good. Contains plotting and optimal path generation capabilities.

## Requirements

Use this file to install:  [Project.toml](Project.toml)
- Interpolations = "0.15.1"
- LinearAlgebra = "1.11.0"
- Plots = "1.40.8"
- Roots = "2.2.1"

## Usage
```julia
include("vfi.jl") # must be placed in working directory
using .vfi

# Define production function 
function F(m :: Model; k)  
    return m.A * k.^m.α
end

# Define consumption function 
function c(m :: Model; k, kPrime)  
    return m.f_fn(m; k = k) .+ (1 - m.δ) .* k .- kPrime
end

# Define CRRA utility function 
function u(m :: Model; c)  
    # Set c to m.lb_c whenever c <= 0
    c = ifelse.(c .<= 0, m.lb_c, c)
    return m.σ == 1 ? log.(c) : ((c .^ (1.0 - m.σ)) .- 1.0) ./ (1.0 - m.σ) 
end

# Define k' constraint
function budget_kprime(m :: Model; c, y, k)  
    return (y .+ (1 - m.δ) .* k) .- c
end

# Define c constraint
function budget_c(m :: Model; kPrime, y, k)  
    return (y .+ (1 - m.δ) .* k) .- kPrime
end

# Set exogenous params
β = 0.96
δ = 0.1
σ = 1.0
A = 1.0
α = 0.33

# Bounds
n = 100
lb_k = 0.05
ub_k = 5.00
lb_kprime = 0.0
v0 = zeros(n)
lb_c = 1e-10
lb_p = -Inf
tol = 1e-6
imax = 1000

# Initialize
results = []

m = vfi.Model(
        β, δ, σ, A, α,
        n, lb_k, ub_k, lb_kprime, v0, lb_c, lb_p, tol, imax,
        F, c, u, budget_kprime, budget_c
    )
    
output = vfi.solver_norisk(m)

push!(results, output)
```

### Plotting the steady state
```julia
vfi.plotter_ss(;df = results, nrow = length(results), ncol = 1, addk = true, var = :n)
```
[plot_ss](example/plot_ss.svg)

### Plotting capital evolution for a given k0 and T
```julia
# Initialize
k0 = 0.1
T = 50
plot_k_path(;k0 = k0, kprime = results[1].kprime, kgrid = results[1].kgrid, T = T)
```
[plot_k_path](example/plot_k_path.svg)

## References
[1]  Permanent Income Model 4: Value Function Iteration https://lhendricks.org/econ890/julia/pih/pih4.html 

[2]  PihVfi890 https://github.com/hendri54/PihVfi890 