# simulate.jl
# -----------------------------------------------------------------------------------------------
# To simulate the economy using the solved model.
# -----------------------------------------------------------------------------------------------

using Distributions

function simulate_model(T_sim, θ::ModelParameters, gp::GridParameters, Dgrid, Y_grid, R_grid, Pi, gD)
    D_sim = zeros(T_sim)
    Y_sim = zeros(T_sim)
    R_sim = zeros(T_sim)
    cT_sim = zeros(T_sim)
    
    D_sim[1] = 0.0 
    state_idx = Int(floor(length(Y_grid)/2)) # to start near steady state
    
    for t in 1:(T_sim-1)
        Y_sim[t] = Y_grid[state_idx]
        R_sim[t] = R_grid[state_idx]
        
        policy_t = gD[:, state_idx]

        D_sim[t+1] = interp1D(Dgrid, policy_t, D_sim[t], gp)
        
        cT_sim[t] = Y_sim[t] - D_sim[t] + D_sim[t+1] / (1.0 + R_sim[t])
        
        p_dist = Pi[state_idx, :]
        state_idx = rand(Categorical(p_dist))
    end
    
    Y_sim[T_sim] = Y_grid[state_idx]
    R_sim[T_sim] = R_grid[state_idx]
    cT_sim[T_sim] = Y_sim[T_sim] - D_sim[T_sim] + D_sim[T_sim] / (1.0 + R_sim[T_sim])
    
    return cT_sim, D_sim, Y_sim, R_sim
end