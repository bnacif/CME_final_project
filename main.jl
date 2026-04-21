# main.jl
# ---------------------------------------------------------------------------
# Computational Methods in Economics - Final Project
# Breno Avila
# April 2026
# ---------------------------------------------------------------------------
# Note: code structure inspired by Professor Lucas Finamor's course codes.
# ---------------------------------------------------------------------------

#=
This codebase solves and simulates an infinite-horizon small open economy 
model based on Schmitt-Grohé and Uribe (2016), extended to include a dual 
labor market with formal and informal sectors. 

The economy is subject to stochastic shocks to both tradable income and 
interest rate.

Version: Optimal Exchange Rate Policy
- The central planner chooses the optimal path of tradable consumption and 
external debt to maximize households expected lifetime utility.
- The solution is approximated using Value Function Iteration (VFI).
=#

#=
ATTENTION!
If this is the first time running the code, run "init_packages.jl" to create the local environment.
This guarantees the use of the same package versions.
=# 

# ---------------------------------------------------------------------------

# Load packages
using Pkg
Pkg.activate(".")    # use the exact versions specified in my Manifest.toml

using Plots
using JLD2
using BenchmarkTools
using LaTeXStrings

include("parameters.jl")
include("functions.jl")
include("procedures.jl")
include("solveModel.jl")
include("simulate.jl")

# Initialize parameters
θ = ModelParameters()
gp = GridParameters() 

Dgrid = get_debt_grid(θ, gp)
Y_grid, R_grid, Pi = get_var1_process(θ, gp)

# Get static optimal non-tradable output
yN_star = get_optimal_yN(θ)
println("Optimal Nontradable Output Level: ", yN_star)
println("State Space Size (Income x Interest Rate): ", gp.NS)

# Solve the Model
# I used Brent for efficiency
V, gD = @btime solveModel($θ, $gp, $Dgrid, $Y_grid, $R_grid, $Pi, $yN_star, opt_algo=Brent())

# Save Data
@save "optimal_policy_var1.jld2" V gD Dgrid Y_grid R_grid θ gp

# Simulate
T_sim = 100
cT_sim, D_sim, Y_sim, R_sim = simulate_model(T_sim, θ, gp, Dgrid, Y_grid, R_grid, Pi, gD)

# Plot
p1 = plot(1:T_sim, cT_sim, label=L"Tradable Consumption ($c^T$)", lw=2)
p2 = plot(1:T_sim, D_sim, label=L"Debt ($d$)", color=:red, lw=2)
p3 = plot(1:T_sim, Y_sim, label=L"Tradable Endowment ($y^T$)", color=:green, lw=2)
p4 = plot(1:T_sim, R_sim .* 100, label=L"Interest Rate ($r$, %)", color=:purple, lw=2)

final_plot = plot(p1, p2, p3, p4, layout=(4,1), size=(800,800))

display(final_plot)
#savefig(final_plot, "simulation_results.png") # save plot to use in .qmd