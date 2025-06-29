using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

include("src/bridge_simulation.jl")

# Beam and material parameters
L = 300.0               # Beam length (m)
n_elem = 25           # Number of finite elements
n_node = n_elem + 1   # Number of nodes
ρ = 7800.0            # Density (kg/m^3)
A = 4.0              # Cross-section area (m^2)
I = 3.0              # Moment of inertia (m^4)
E0 = 207e9             # Base Young's modulus (Pa)
α = -1e5              # E-temperature slope (Pa/K)
cutoff_freq = 25.0  # Cutoff frequency for modes (Hz)

bc = BridgeBC([  # Node 1: both translational and rotational DOFs fixed
    [1, "trans"],
    [n_node, "y"]
])

bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 206e9; 70 150e9], cutoff_freq)

@info "Number of modes: $n_modes"

# Support element parameters with temperature dependence
E_bridge = [
    -10.0  250e9;   # E at -10°C
     20.0  207e9;   # E at 20°C
     50.0  150e9    # E at 50°C (thermal expansion reduces stiffness)
]

E_support = E_bridge

A_support = 0.01        # Cross-section area (m^2)
I_support = 0.001       # Moment of inertia (m^4)
L_support = 50.0       # Length of support element (m)

se = [SupportElement(
        n_elem ÷ 3,          # Connect to middle of bridge
        [1, 2, 3],               # Connect x,y,ϕ DOFs
        -150.0,                   # angle (degrees)
        5,                       # 5 elements in support
        A_support,               # Cross-sectional area
        I_support,               # Moment of inertia
        E_support,          # Temperature-dependent Young's modulus
        L_support,               # Length
        BCTypes["trans"]           # Fix all DOFs at bottom
    ),
    SupportElement(
        n_elem - (n_elem ÷ 3),          # Connect to middle of bridge
        [1, 2, 3],               # Connect x,y,ϕ DOFs
        -30.0,                   # angle (degrees)
        5,                       # 5 elements in support
        A_support,               # Cross-sectional area
        I_support,               # Moment of inertia
        E_support,          # Temperature-dependent Young's modulus
        L_support,               # Length
        BCTypes["trans"]           # Fix all DOFs at bottom
    )]

# Example workflow with supports:
supports = se  # Your support element
# supports = SupportElement[]

Ts = [-10.0, 20.0, 50.0]  # Temperatures in degrees Celsius

# Create comprehensive simulation options
sim_opts = SimulationOptions(
    bo, supports, collect(Ts), damping_ratio=0.02
)
save_simulation_options(sim_opts, "data/bridge_simulation_config.json")

# supports = [se]  # Add your support element
M, K, λs, vectors, vectors_unnormalized = assemble_and_decompose(bo, collect(Ts); supports=supports)

# Solve dynamics (same ODE, but with expanded system)
n_modes = size(λs, 1)
total_dofs = size(vectors, 1)

# Create interpolators (same as before)
ω_interp = interpolate((1:size(λs,1), Ts), λs, Gridded(Linear()))
Φ_interp = interpolate((1:size(vectors,1), 1:size(vectors,2), Ts), vectors, Gridded(Linear()))

ω_T = t -> ω_interp(1:n_modes, t)
Φ_T = t -> Φ_interp(1:total_dofs, 1:n_modes, t)

# Your load_vector function needs to account for expanded DOF numbering
load_vector = (t, dof) -> begin 
    f = zeros(length(dof))
    # Apply load to expanded system...
    return f
end

force_node = collect(2:3:3*bo.n_elem)  # Node to apply force on (y-displacement)

tspan = (0.0, 300.0)

load_vector = (t, dof) -> begin 
    f = zeros(length(dof))
    
    # Linear frequency sweep from f1 to f2
    f1 = 0.1  # Starting frequency (Hz)
    f2 = 10.0 # Ending frequency (Hz)
    
    # Instantaneous frequency: f(t) = f1 + (f2-f1) * t/T
    freq_t = f1 + (f2 - f1) * (t / tspan[2])
    
    # Phase accumulation for chirp: φ(t) = 2π ∫₀ᵗ f(τ) dτ
    phase = 2π * (f1 * t + (f2 - f1) * t^2 / (2 * tspan[2]))
    
    y_dof = force_node
    f[y_dof] .= 1000.0 * sin(phase)
    
    return f
end

damping_ratio = 0.0000

u0 = zeros(2 * n_modes)

T_func = (t) -> Ts[end] - (Ts[end] - Ts[1]) * (t / tspan[2])  # Linear temperature change from Ts[1] to Ts[end]

prob = ODEProblem(beam_modal_ode!, u0, tspan,
                    (; 
                        n_modes=n_modes,
                        n_dofs=total_dofs,  # FIXED: Use total DOFs including supports
                        T_func = T_func,
                        ω_interp = ω_T,
                        ζ = fill(damping_ratio, n_modes),  # Constant damping
                        Φ_interp = Φ_T,
                        load_vector = load_vector,
                    ))

@info "Solving dynamic response..."
@time sol = solve(prob, BS3(), saveat=0.01)

q = reduce(hcat, sol.u)

u, du = reconstruct_physical(sim_opts, q, Φ_T, T_func, sol.t)
# # 1. Plot structure only
# plot_bridge_with_supports(bo, supports)

time_subsample = sol.t[1:2:end]
u_subsample = u[:, 1:2:end]
anim_fast = animate_dynamic_response(bo, supports, u_subsample, time_subsample,
                                   scale_factor=100.0,
                                   n_frames=500,
                                   fps=24,
                                   filename="bridge_fast_dynamics.gif")

# support_dof_maps, total_dofs = create_support_dof_mapping(bo, supports)

# 2. Plot specific mode shape
mode_num = 3
mode_shape = Φ_T(50.0)[:,mode_num]  # First mode at first temperature

plot_mode_shape(bo, supports, mode_shape, mode_num, scale_factor=10000.0)

# # 3. Animate a mode
# mode_num = 1
# mode_shape = vectors_unnormalized[:, mode_num, 1]  # Second mode
# anim = animate_mode(bo, supports, mode_shape, mode_num, scale_factor=10.0, fsize=(800, 600))
# gif(anim, "mode2_animation.gif", fps=15)