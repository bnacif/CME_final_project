# parameters.jl 
# -----------------------------------------------------------------------------------------------
# To define parameters used in the model.
# -----------------------------------------------------------------------------------------------

using Parameters

# Parameters
# -----------------------

# Most parameters follow those calibrated and estimated using Argentine data in Schmitt-Grohé and Uribe (2016).
# Specifically, the paper sets the labor share in the nontraded sector α = 0.75. 
# Since this new version considers a dual labor market, I split this total share into formal (α_F = 0.50) and informal (α_I = 0.25$) 
# components to preserve the aggregate labor returns to scale.

# Note: @with_kw is a macro from Parameters package. It allows a keyword constructor with default values (that can be changed individually).

@with_kw struct ModelParameters
    # Preferences
    β::Float64 = 0.9375       # Discount factor
    σ::Float64 = 5.0          # Inverse of intertemporal elasticity
    ξ::Float64 = 0.44         # Elasticity of sub between tradable/nontradable
    a::Float64 = 0.26         # Share of tradables in CES
    
    # Technology & Labor
    h_bar::Float64 = 1.0      # Aggregate labor endowment
    α_F::Float64 = 0.50       # Formal labor share
    α_I::Float64 = 0.25       # Informal labor share
    
    # Debt Constraints
    d_max::Float64 = 8.34     # Natural debt limit (borrowing constraint)
    d_min::Float64 = -5.0     # Maximum savings limit
    
    # SGU (2016) VAR(1) Estimates for Argentina
    # State vector is [ln(y_t), ln((1+r_t)/(1+r_ss))]
    A::Matrix{Float64} = [0.7901 -1.3570; -0.0104 0.8638]
    Ω::Matrix{Float64} = [0.0012346 -0.0000776; -0.0000776 0.0000401]
    r_ss::Float64 = 0.0316  # Steady state interest rate
    
    # Minimum consumption
    minCons::Float64 = 1e-4
end


# Grid
# -----------------------

# Changing the number of points in the each grid significantly changes computational time.
# Below, I repeat the values in SGU (2016).
# Still, a corser discretization yileds a similar policy function

@with_kw struct GridParameters
    ND::Int = 50             # Grid points for debt (leads the increase in computational time)
    NY::Int = 21              # Grid points for income
    NR::Int = 11              # Grid points for interest rate
    NS::Int = NY * NR         # Total states (231)
    tol::Float64 = 1e-6       # VFI Tolerance
    max_iter::Int = 2000      # VFI Max Iterations

    # computational methods
    grid_spacing::String = "equalsteps"   # "equalsteps" or "logsteps"
    interp_method::String = "linear"      # "linear" or "pchip"
end
