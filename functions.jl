# functions.jl
# -----------------------------------------------------------------------------------------------
# To define main functions used in the project.
# -----------------------------------------------------------------------------------------------

using LinearAlgebra
using Distributions

# Static optimal nontradable output (first-best allocation)
function get_optimal_yN(θ::ModelParameters)
    share_F = θ.α_F / (θ.α_F + θ.α_I)
    share_I = θ.α_I / (θ.α_F + θ.α_I)
    
    h_F = share_F * θ.h_bar
    h_I = share_I * θ.h_bar
    
    yN_star = (h_F^θ.α_F) * (h_I^θ.α_I)
    return yN_star
end

# CES Aggregator
function aggregator(cT, cN, θ::ModelParameters)
    term1 = θ.a * (cT)^(1.0 - 1.0/θ.ξ)
    term2 = (1.0 - θ.a) * (cN)^(1.0 - 1.0/θ.ξ)
    return (term1 + term2)^(θ.ξ / (θ.ξ - 1.0))
end

# Utility Function
function utility(cT, cN, θ::ModelParameters)
    if cT < θ.minCons
        return -1e10 # Heavy penalty for violating constraints
    end
    c = aggregator(cT, cN, θ)
    return (c^(1.0 - θ.σ) - 1.0) / (1.0 - θ.σ)
end

# Discretization of y^T and r processes
function get_var1_process(θ::ModelParameters, gp::GridParameters; m=3.0)
    # Compute unconditional variance via discrete Lyapunov approx
    V = copy(θ.Ω)
    for i in 1:500
        V = θ.A * V * transpose(θ.A) + θ.Ω
    end
    
    std_y = sqrt(V[1,1])
    std_r = sqrt(V[2,2])
    
    # Create grids
    y_log_grid = collect(range(-m*std_y, m*std_y, length=gp.NY))
    r_log_grid = collect(range(-m*std_r, m*std_r, length=gp.NR))
    
    # Build combined state space
    states = zeros(gp.NS, 2)
    idx = 1
    for ir in 1:gp.NR
        for iy in 1:gp.NY
            states[idx, 1] = y_log_grid[iy]
            states[idx, 2] = r_log_grid[ir]
            idx += 1
        end
    end
    
    # Build Transition Matrix (Pi)
    Pi = zeros(gp.NS, gp.NS)
    dist = MvNormal(zeros(2), θ.Ω)
    
    for i in 1:gp.NS
        mu = θ.A * states[i, :] # Expected next state
        for j in 1:gp.NS
            # PDF midpoint approximation for transition probability
            Pi[i, j] = pdf(dist, states[j, :] .- mu)
        end
        # Normalize row to sum to 1
        Pi[i, :] ./= sum(Pi[i, :])
    end
    
    # Convert states back to levels
    Y_grid = exp.(states[:, 1])
    R_grid = (1.0 + θ.r_ss) .* exp.(states[:, 2]) .- 1.0
    
    # Return vectors of length NS, mapped to the NSxNS transition matrix
    return Y_grid, R_grid, Pi
end