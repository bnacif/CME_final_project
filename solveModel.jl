# solveModel.jl
# -----------------------------------------------------------------------------------------------
# To define functions that solve the model using the equilibrium equations.
# -----------------------------------------------------------------------------------------------

using Optim
using Base.Threads

function solveModel(θ::ModelParameters, gp::GridParameters, Dgrid, Y_grid, R_grid, Pi, yN_star; opt_algo=Brent())
    # Empty objects to store functions
    V = zeros(gp.ND, gp.NS)
    V_new = zeros(gp.ND, gp.NS)
    gD = zeros(gp.ND, gp.NS)
    
    err = 1.0
    iter = 0
    
    while err > gp.tol && iter < gp.max_iter
        iter += 1
        
        # VFI
        # Multi-threading (parallelization of the outer loop)
        Threads.@threads for is in 1:gp.NS
            Y_now = Y_grid[is]
            R_now = R_grid[is] 
            
            # Note: Because EV_next is allocated inside the threaded loop, each thread gets its own local copy
            EV_next = zeros(gp.ND)
            for id_next in 1:gp.ND
                EV_next[id_next] = sum(Pi[is, is_next] * V[id_next, is_next] for is_next in 1:gp.NS)
            end
            
            for id in 1:gp.ND
                D_now = Dgrid[id]
                
                function objFunc(D_next)
                    cT = Y_now - D_now + D_next / (1.0 + R_now)
                    if cT < θ.minCons
                        return 1e10    # penalty to guarantee positive consumption
                    end
                    # Interpolate the expected value function
                    # using the method in parameters.jl
                    val_next = interp1D(Dgrid, EV_next, D_next, gp)
                    return -(utility(cT, yN_star, θ) + θ.β * val_next)  # - for minimization
                end
                
                res = optimize(objFunc, θ.d_min, θ.d_max, opt_algo)
                
                gD[id, is] = res.minimizer
                V_new[id, is] = -res.minimum
            end
        end
        
        err = maximum(abs.(V_new .- V))
        V .= V_new
        
        if iter % 50 == 0
            println("Iteration: $iter, Error: $err")
        end
    end
    
    println("VFI Converged in $iter iterations.")
    return V, gD
end