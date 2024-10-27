"""
Generic solver for the optimal neoclassical growth model, using value function iteration (VFI).
"""
module vfi

using Interpolations, Roots, LinearAlgebra, Plots

export Model, solver_norisk, plotter_ss, gen_k_path, plot_k_path

struct Model
    # parameters
    β::Float64         # Discount factor    
    δ::Float64          # Depreciation rate
    σ::Float64          # Risk aversion  
    A::Float64          # Technology parameter
    α::Float64        # Capital share
    # Bounds
    n_k::Int64      # Number of grid points for capital
    lb_k::Float64   # Lower bound of capital grid
    ub_k::Float64   # Upper bound of capital grid
    lb_kprime::Float64           # minimum investment (if not endogenous)
    v_0:: Vector{Float64}   # Initial value function
    lb_c::Float64   # Lower bound of consumption grid
    lb_p::Float64   # Lower bound of return space
    tol::Float64            # Convergence error tolerance
    imax::Int64             # Maximum number of iterations
    # functions
    f_fn::Function             # Production function
    c_fn::Function                # Consumption function
    u_fn::Function                # Utility function
    b_kprime_fn::Function                # Budget constraint function for
    b_c_fn::Function                # Budget constraint function for consumption
end

## ----------  Helpers
# Define function to generate discrete state space
get_k_grid(m :: Model) = collect(LinRange(m.lb_k, m.ub_k, m.n_k));

## ---------- Constraints Handling
# Validation check: Max kPrime consistent with c > c_min
function kprime_max(m :: Model; k_grid)
    m.b_kprime_fn(m; c = m.lb_c, y = m.f_fn(m, k = k_grid), k = k_grid)
end;

function kprime_min(m :: Model; k_grid)
    m.b_kprime_fn(m; c = k_grid, y = 0.0, k = k_grid)
end;

## ---------- Value Function Iteration (VFI) Solver

"""
solve_vfi: Solve for an optimal policy using value function iteration (VFI).
"""
function solver_norisk(m::Model)
    # Initialize
    N = 0 # Terminal step
    T = m.n_k # Grid size
    kp = Vector{Any}(undef, T);
    c_path = Vector{Any}(undef, T);
    β = m.β
    σ = m.σ
    tol = m.tol
    imax = m.imax
    v = m.v_0
    v_new = deepcopy(v)
    error = Inf

    # Generate grid
    k = get_k_grid(m)

    # Generate indices for undefined T x 1 array
    idxs = Array{CartesianIndex{2}, 1}(undef, T)

    # Consumption 
    cgrid = m.c_fn(m; k = k, kPrime = k') # T x T matrix of consumption

    # Apply boundary conditions for valid k'
    kp_min = m.lb_kprime
    kp_max = m.b_kprime_fn(m; c = cgrid, y = m.f_fn(m; k = k), k = k)

    # Mask invalid k' choices
    valid_mask = (k .>= kp_min) .& (k .<= kp_max)
    invalid_mask = .!valid_mask .* .!valid_mask' 

    # Compute utility for valid consumption only
    ugrid = m.u_fn(m; c = cgrid)  # T x T matrix of utility
    ugrid[invalid_mask] .= m.lb_p  # Invalid choices set to lower bound of return space

    # Backward induction loop
    while N <= imax && error > tol
        N += 1 
        v_interpol = LinearInterpolation(k, v) # allowing kprime to not be on grid

        # Calculate the utility at each grid point
        obj = ugrid .+ β .* v_interpol(k)'


        # Optimize k' to max u
        if σ > 1.0 # CRRA not bounded from below if σ > 1
            (v_new, idxs) = findmin(-1.0.*obj, dims=2) # flip sign to find max
            v_new = -v_new
        else
            (v_new, idxs) = findmax(obj, dims=2)
        end

        error = maximum(abs.(v_new - v))
        # print("\nIteration: $(N), Error: $(error)")

        # Bellman update for value function
        v .= v_new
    end

    # Extract law of motion for k
    kprimeidx = getindex.(idxs, 2)  
    kprime = dropdims(k[kprimeidx], dims = 2)

    # Compute the optimal consumption policy given k and k'
    c_path = m.b_c_fn(m; kPrime = kprime, y = m.f_fn(m, k = k), k = k);

    print("\nConverged after $(N) iterations.")

    return (n = T, N = N, kgrid = k, v = v, kprime = kprime, kprimeidx = kprimeidx, imax = imax, σ = σ, pol_c = c_path)
end

function plotter_ss(;df, nrow = 1, ncol = 1, addk = false, var::Symbol)
    #Initialize
    p = []
    all_x = Any[]  
    all_y = Any[]  
    all_lb = String[]

    for i in eachindex(df)
        plot_df = df[i]
        x = plot_df.kgrid
        val = getfield(plot_df, var)
        title_str = "$(var) = $(val)"
        
        if addk == true
            y = hcat(plot_df.v, plot_df.kprime) 
            lb = ["v(k)" "k'(k)"]
            plot = scatter(
                x, y, ms=2, ma=0.5, msw =0.2, title = title_str,
                label= lb , xlabel="k", legend=false
            )
            push!(p, plot)
        else
            push!(all_x, x)  
            push!(all_y, plot_df.v)  
            push!(all_lb, title_str)
        end
    end
    
    if addk == true
        p[nrow] = plot!(p[nrow], legend = :bottomright)
        plot(p..., layout=(nrow, ncol), size = (1000, 1000))
    else
        colors = palette(:tab10, nrow) 
        markers = filter((m->begin
                m in Plots.supported_markers()
            end), Plots._shape_keys)
        plot()
        for i in 1:length(all_x)
            scatter!(all_x[i], all_y[i], m = markers[i], ms=5, ma=0.2, msw =0.2, msa =0.9, label=all_lb[i], color=colors[i], xlabel="k", legend=true, size = (1000, 500))
        end
        plot!(legend=:outerbottom, legendcolumns=nrow)
        current()
    end
end


function gen_k_path(; k0, kprime, kgrid, T)
    kprime_interp = LinearInterpolation(kgrid, kprime, extrapolation_bc=Line()) 
    capital_path = fill(k0, T)

    for t in 2:T
        capital_path[t] = kprime_interp(capital_path[t-1])
    end

    return capital_path
end

function plot_k_path(; k0, kprime, kgrid, T)
    k_path = gen_k_path(; k0 = k0, kprime = kprime, kgrid = kgrid, T = T)

    plot()
    plot!(0:T-1, k_path, xlabel="Time (t)", ylabel="Capital (kₜ)", title="Capital Path over Time", legend = false)
    annotate!([(T, k_path[T], ("k_T = $(round(k_path[T]; digits=2))", 8, :red, :right))])
    current()
end

end  # module vfi
