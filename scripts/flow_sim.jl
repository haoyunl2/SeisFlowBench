## ============================================================================
## Immiscible Two-Phase Flow Simulations: CO₂ and Water
##
## This script simulates immiscible two-phase flow in a porous medium using
## the JutulDarcyRules.jl framework. The simulation models CO₂ injection and
## its displacement of water in the subsurface, tracking pressure evolution
## and saturation changes.
##
## Dependencies:
## - JutulDarcyRules.jl: Reservoir simulation package
## - Jutul.jl: Computational framework for PDE-constrained optimization
## - JutulDarcy.jl: Darcy flow simulation module
## - DrWatson.jl: Experiment management
## - JLD2.jl: Data storage
## - PyPlot.jl: Visualization
##
## Output:
## - Saturation (`S_array`)
## - Pressure (`P_array`)
## - Pressure difference (`P_diff_array`)
## - Bottom-hole pressure (`BHP_array`)
## - Plots of permeability, saturation, and pressure difference
##
## Simulation results are stored in JLD2 format and can be converted to HDF5.
##
## Author: [Haoyun Li and Abhinav Gahlot]
## Date: [Feb 25th, 2025]
## ============================================================================

using Pkg
Pkg.activate(".")  # Activate the project environment

# Uncomment this line when running for the first time to install dependencies
# Pkg.instantiate()

using DrWatson, JutulDarcyRules, LinearAlgebra, JLD2, PyPlot, PyCall, SlimPlotting
@pyimport cmasher  # Colormap utilities

# -----------------------------------------------------------------------------
# Domain Parameters
# -----------------------------------------------------------------------------
n = (128, 1, 128)  # Grid dimensions
d = (6.25 * 2, 100.0, 6.25 * 2)  # Grid spacing (dx, dy, dz)
h = 0.0  # Initial depth reference
ϕ = 0.25  # Porosity

# Compute initial pressure field (hydrostatic pressure)
P0 = (repeat(collect(1:128), 1, 128) * d[3] .+ h) * JutulDarcyRules.ρH2O * 9.807

# -----------------------------------------------------------------------------
# Simulation Settings
# -----------------------------------------------------------------------------
num_sample = 10  # Number of permeability samples to simulate
perturb = true  # Use perturbed permeability fields

# Load permeability dataset
perm_file = perturb ? "K_array_perturb.jld2" : "K_array_wise.jld2"
JLD2.@load datadir("perm", perm_file) K_array

# -----------------------------------------------------------------------------
# Simulation Function
# -----------------------------------------------------------------------------
"""
    simulation(irate, K, tstep)

Runs the two-phase flow simulation with CO₂ injection.

# Arguments:
- `irate`: Injection rate (CO₂)
- `K`: Permeability field
- `tstep`: Time-stepping array

# Returns:
- `S_array`: Saturation over time
- `P_array`: Pressure over time
- `P_diff_array`: Pressure difference over time
- `BHP_array`: Bottom-hole pressure at the injector
"""
function simulation(irate, K, tstep)
    # Determine injection location based on highest permeability zone
    inj_GT_2D = 101 + argmax(@view K[63, 101:110]) - 1
    inj_loc_grid = (63, 1, inj_GT_2D)
    inj_loc = inj_loc_grid .* d

    num_t = length(tstep)

    # Set up the reservoir model
    model = jutulModel(n, d, ϕ, K1to3(K; kvoverkh=0.36); h=h)
    f = jutulVWell(irate, [(inj_loc[1], inj_loc[2])]; startz=[inj_loc[3]], endz=[inj_loc[3] + 6 * d[3]])
    S = jutulModeling(model, tstep)
    Trans = KtoTrans(CartesianMesh(model), K1to3(K; kvoverkh=0.36))

    # Run the simulation
    @time states = S(log.(Trans), f)

    # Extract results
    S_array = [reshape(states.states[i][1:n[1]*n[3]], n[1], n[end]) for i in 1:num_t]
    P_array = [reshape(states.states[i][n[1]*n[3]+1:end], n[1], n[end]) for i in 1:num_t]
    P_diff_array = [@views reshape(states.states[i][n[1]*n[3]+1:end], n[1], n[end]) - P0' for i in 1:num_t]
    BHP_array = [states.states[i].state[:Injector][:Pressure] for i in 1:num_t]

    GC.gc()
    return S_array, P_array, P_diff_array, BHP_array
end

# -----------------------------------------------------------------------------
# Visualization Functions
# -----------------------------------------------------------------------------

"""
    plot_permeability(K, sample, path)

Generates and saves a permeability plot.
"""
function plot_permeability(K, sample, path)
    fig = figure(figsize=(5, 5))
    imshow(K' / md, extent=(0, (n[1] - 1) * d[1], h + (n[3] - 1) * d[3], h))
    clb_K = colorbar()
    clb_K[:ax][:set_title]("Md", fontsize=12)
    xlabel("X [m]", fontsize=12)
    ylabel("Depth [m]", fontsize=12)
    title("Permeability", fontsize=15)
    plt.tight_layout()
    
    safesave(joinpath(path, savename(@strdict(sample); digits=6) * "_perm.png"), fig)
    close(fig)
end

"""
    plot_saturation_pressure(S_array, P_diff_array, sample, path, tstep)

Generates and saves saturation and pressure difference plots over time.
"""
function plot_saturation_pressure(S_array, P_diff_array, sample, path, tstep)
    st, tl = Int(length(tstep) / 10), length(tstep)  # Step interval for plotting

    for step in st:st:tl
        fig = figure(figsize=(5, 10))

        # Saturation plot
        subplot(2, 1, 1)
        imshow(S_array[step]', vmin=0, vmax=1)
        colorbar()
        title("Saturation", fontsize=15)
        xlabel("X[m]")
        ylabel("Depth[m]")

        # Pressure difference plot
        subplot(2, 1, 2)
        imshow(P_diff_array[step]')
        clb_P = colorbar()
        clb_P[:ax][:set_title]("MPa", fontsize=12)
        title("Pressure Difference", fontsize=15)
        xlabel("X[m]")
        ylabel("Depth[m]")

        safesave(joinpath(path, savename(@strdict(sample, step); digits=6) * "_sat&presdiff.png"), fig)
        close(fig)
    end
end

# -----------------------------------------------------------------------------
# Main Simulation Loop
# -----------------------------------------------------------------------------
tstep = 8 * ones(6 * 2 * 10)  # Time steps

sim_name = perturb ? "2D_perm_perturb" : "2D_perm_wise"
exp_name = "flow"  # Define a variable first
plot_path = plotsdir(sim_name, savename(@strdict(exp_name); digits=6))


for sample in 1:num_sample
    irate = 0.05
    K = K_array[sample, :, :] * JutulDarcyRules.md

    # Plot permeability field
    plot_permeability(K, sample, plot_path)

    # Run the simulation
    S_array, P_array, P_diff_array, BHP_array = simulation(irate, K, tstep)

    # Plot saturation and pressure difference
    plot_saturation_pressure(S_array, P_diff_array, sample, plot_path, tstep)

    # Save results
    @tagsave(datadir(sim_name, savename(@strdict(sample), "jld2"; digits=6)),
        Dict("S_array" => S_array, "P_array" => P_array, "P_diff_array" => P_diff_array, "BHP_array" => BHP_array);
        safe=true)
end

println("✅ Simulation and data saving complete!")
