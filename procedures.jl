# procedures.jl
# -----------------------------------------------------------------------------------------------
# To define general routines used in the project.
# -----------------------------------------------------------------------------------------------

using Interpolations

# Debt grid
# -----------------------

# In SGU(2016), all grids are equally spaced.

function get_debt_grid(θ::ModelParameters, gp::GridParameters)    
    # initialize the grid
    Dgrid = zeros(gp.ND)
    Dlog  = zeros(gp.ND)

    # set minimum and maximum values
    Dgrid[1]  = θ.d_min
    Dgrid[gp.ND] = θ.d_max

    # compute span 
    Span = θ.d_max - θ.d_min

    # Replace aux grid to be in logs
    Dlog[1] = 0.0 
    Dlog[gp.ND] = Span
    
    # Build the grid based on the chosen spacing method
    # In SGU(2016), all grids are equally spaced, but we allow for logsteps 
    # to concentrate points near the borrowing constraint if needed.
    for j in 2:gp.ND
        if (gp.grid_spacing == "equalsteps")
            Dlog[j] = Dlog[j-1] + Span/(gp.ND-1)
        elseif (gp.grid_spacing == "logsteps")
            Dlog[j] = Dlog[j-1] + log(Span+1.0)/(gp.ND-1)
        else 
            throw("grid_spacing should be either equalsteps or logsteps")
        end
    end

    # Transform the grid back to levels
    for j in 1:gp.ND
        if (gp.grid_spacing == "equalsteps")
            Dgrid[j] = Dlog[j] + θ.d_min
        elseif (gp.grid_spacing == "logsteps")
            Dgrid[j] = θ.d_max - (exp(Dlog[gp.ND - j + 1]) - 1.0) # dense points at the borrowing constraint (d_max)
        end
    end

    return Dgrid
end


# Interpolation 
# -----------------------

# In SGU(2016), there is no interpolation, but the continuous choice
# in this version requires this addition.

function interp1D(Xgrid, Ygrid, x, gp::GridParameters)

    if (gp.interp_method == "linear")
        itp = linear_interpolation(Xgrid, Ygrid)
        return itp(x)

    elseif (gp.interp_method == "pchip")
        itp = interpolate(Xgrid, Ygrid, FritschButlandMonotonicInterpolation())
        return itp(x)

    else
        throw("interp_method should be either linear or pchip")
    end 
end
