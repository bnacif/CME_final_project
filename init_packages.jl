# init_packages.jl
# -----------------------------------------------------------------------------------------------
# To initialize all packages required to run the project
# Run it only once, before main.jl 
# -----------------------------------------------------------------------------------------------

using Pkg
Pkg.activate(".")  # Activates the project environment (Project.toml) in the current directory

# Add required packages in the project environment
Pkg.add([
    "Optim", 
    "Plots", 
    "Interpolations",
    "Parameters", 
    "QuantEcon", 
    "Distributions",
    "JLD2",
    "BenchmarkTools",
    "LaTeXStrings"
])

Pkg.instantiate()  # Reads Manifest.toml and fetches the exact versions
Pkg.status()       # Check package status

# In main.jl and other files, the code calls each package in its correct version
