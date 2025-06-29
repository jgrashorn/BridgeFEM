using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

include("src/bridge_simulation.jl")

# Beam and material parameters
L = 300.0               # Beam length (m)
n_elem = 20           # Number of finite elements
n_node = n_elem + 1   # Number of nodes
ρ = 7800.0            # Density (kg/m^3)
A = 4.0              # Cross-section area (m^2)
I = 3.0              # Moment of inertia (m^4)
E0 = 207e9             # Base Young's modulus (Pa)
α = -1e5              # E-temperature slope (Pa/K)
cutoff_freq = 100.0  # Cutoff frequency for modes (Hz)

bc = BridgeBC([  # Node 1: both translational and rotational DOFs fixed
    [1, "trans"],
    [n_node, "y"]
])

bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 206e9; 70 150e9], cutoff_freq)

save_options_to_json(bo, T=collect(0.0:10.0:40.0))

nTs = 5
Ts = range(10.0, stop=30.0, length=nTs)  # Temperatures to sample
M, K, λs, vectors = assemble_and_decompose(bo, collect(Ts))

n_modes = size(λs, 1)

@info "Number of modes: $n_modes"

ω_interp = interpolate((1:n_modes, Ts), λs, Gridded(Linear()))
Φ_interp = interpolate((1:bo.n_dofs, 1:n_modes, Ts), vectors, Gridded(Linear()))

# Time span
tspan = (0.0, 10.0)

T_func = t -> Ts[end] + (Ts[1] - Ts[end]) * t / tspan[2]
load_vector = (t, dof) -> begin 
    f = zeros(length(dof))
    # Apply load to y-displacement (DOF 2) of middle node (node n_elem÷2 + 1)
    middle_node = n_elem ÷ 2 + 1
    y_dof = 3 * (middle_node - 1) + 2  # y-displacement DOF
    f[y_dof] = 1e6 * sin(2π * (t/30.)*10 * t / 10.0)
    return f
end

# Initial conditions: zero displacement and velocity
u0 = zeros(2 * n_modes)
# ODE Problem setup
@info "Setting up ODE problem for beam dynamics"
prob = ODEProblem(beam_modal_ode!, u0, tspan,
                (; 
                    n_modes=n_modes,
                    n_dofs=bo.n_dofs,
                    T_func = T_func,  # Linear interpolation of temperature
                    ω_interp = ω_interp,
                    ζ = zeros(n_modes),  # Damping ratios (could be constant or interpolated)
                    Φ_interp = Φ_interp,
                    load_vector = load_vector,
                ))

# @time sol = solve(prob, BS3(), saveat=0.01)
# q = reduce(hcat,sol.u)

# @time u, du = reconstruct_physical(bo, q, Φ_interp, T_func, sol.t)
# n_dof = bo.n_dofs

# plot_dof = collect(2:3:n_dof) # Middle node for visualization
# plot(sol.t, u[plot_dof, :]', color=:redsblues, xlabel="Time (s)", ylabel="Y-Displacement (m)", title="Beam Dynamics Simulation",
#     legend=:topright)

# Support element parameters
E = 207e9              # Young's modulus (Pa)
A_support = 0.1        # Cross-section area (m^2)
I_support = 0.01        # Moment of inertia (m^4)
L_support = 50.0       # Length of support element (m)

se = [SupportElement(
        n_elem ÷ 2 + 1,          # Connect to middle of bridge
        [1, 2,3],     # Connect x,y,ϕ DOFs
        -90.0,       # angle (degrees)
        5,          # 5 elements in support
        E*A_support, 
        E*I_support,
        L_support,
        BCTypes["all"]   # Fix all DOFs at bottom
    )]

# Example workflow with supports:
supports = se  # Your support element
# supports = SupportElement[]

# supports = [se]  # Add your support element
M, K, λs, vectors, vectors_unnormalized = assemble_and_decompose(bo, collect(Ts); supports=supports)
M, K, λs, vectors, vectors_unnormalized = assemble_and_decompose(so.bridge, so.temperatures; supports=so.supports)
# Create interpolators (same as before)
ω_interp = interpolate((1:size(λs,1), Ts), λs, Gridded(Linear()))
Φ_interp = interpolate((1:size(vectors,1), 1:size(vectors,2), Ts), vectors, Gridded(Linear()))

# Solve dynamics (same ODE, but with expanded system)
n_modes = size(λs, 1)
total_dofs = size(vectors, 1)  # This is now larger than bo.n_dofs

# Your load_vector function needs to account for expanded DOF numbering
load_vector = (t, dof) -> begin 
    f = zeros(length(dof))
    # Apply load to expanded system...
    return f
end

# 1. Plot structure only
plot_bridge_with_supports(bo, supports)

support_dof_maps, total_dofs = create_support_dof_mapping(bo, supports)

# 2. Plot specific mode shape
mode_num = 1
mode_shape = vectors_unnormalized[:, mode_num, 1]  # First mode at first temperature
plot_mode_shape(bo, supports, mode_shape, mode_num, scale_factor=10.0)

# 3. Animate a mode
mode_num = 1
mode_shape = vectors_unnormalized[:, mode_num, 1]  # Second mode
anim = animate_mode(bo, supports, mode_shape, mode_num, scale_factor=10.0, fsize=(800, 600))
gif(anim, "mode2_animation.gif", fps=15)

ω_interp = (p,t) -> λs[:, 1]
Φ_interp = (p,c,t) -> vectors[:,:,1]

force_node = collect(2:3:size(vectors,1))  # Node to apply force on (y-displacement)

load_vector = (t, dof) -> begin 
    f = zeros(length(dof))
    
    # Apply sinusoidal force to y-displacement of specified node
    y_dof = force_node  # y-displacement DOF
    
    # Only apply if this DOF exists in the expanded system
    f[y_dof] .= 10000.0 * sin(2π * 10.0 * t)
    
    return f
end

damping_ratio = 0.0001

tspan = (0.0, 10.0)
u0 = zeros(2 * n_modes)

prob = ODEProblem(beam_modal_ode!, u0, tspan,
                    (; 
                        n_modes=n_modes,
                        n_dofs=total_dofs,  # FIXED: Use total DOFs including supports
                        T_func = T_func,
                        ω_interp = ω_interp,
                        ζ = fill(damping_ratio, n_modes),  # Constant damping
                        Φ_interp = Φ_interp,
                        load_vector = load_vector,
                    ))
    
@info "Solving dynamic response..."
@time sol = solve(prob, BS3(), saveat=0.01)

# Create comprehensive simulation options
sim_opts = SimulationOptions(
    bo, supports, collect(Ts), damping_ratio=0.02
)

save_simulation_options(sim_opts, "data/bridge_simulation_config.json")

so = load_simulation_options("data/bridge_simulation_config.json")

M, K, λs, vectors, vectors_unnormalized = assemble_and_decompose(so.bridge, so.temperatures; supports=so.supports)

q = reduce(hcat, sol.u)