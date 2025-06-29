using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

include("src/bridge_simulation.jl")

# Beam and material parameters
L = 100.0               # Beam length (m)
n_elem = 5           # Number of finite elements
n_node = n_elem + 1   # Number of nodes
ρ = 7800.0            # Density (kg/m^3)
A = 4.0              # Cross-section area (m^2)
I = 3.0              # Moment of inertia (m^4)
E0 = 207e9             # Base Young's modulus (Pa)
α = -1e5              # E-temperature slope (Pa/K)
cutoff_freq = 100.0  # Cutoff frequency for modes (Hz)

bc = BridgeBC([  # Node 1: both translational and rotational DOFs fixed
    [1, "all"],
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

plot(λs')
ylabel!("Frequency (Hz)")
xlabel!("Temperature (°C)")

# Time span
tspan = (0.0, 30.0)

T_func = t -> Ts[end] + (Ts[1] - Ts[end]) * t / tspan[2]
load_vector = (t, dof) -> begin 
    f = zeros(length(dof))
    # Apply load to y-displacement (DOF 2) of middle node (node n_elem÷2 + 1)
    middle_node = n_elem ÷ 2 + 1
    y_dof = 3 * (middle_node - 1) + 2  # y-displacement DOF
    f[y_dof] = 1e6 * sin(2π * 10 * t / 10.0)
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

@time sol = solve(prob, BS3(), saveat=0.01)
q = reduce(hcat,sol.u)

@time u, du = reconstruct_physical(bo, q, Φ_interp, T_func, sol.t)

plot_dof = 3*(n_elem ÷ 2) + 1 # Middle node for visualization
plot(sol.t, u[plot_dof, :], color=:redsblues, xlabel="Time (s)", ylabel="Y-Displacement (m)", title="Beam Dynamics Simulation",
    legend=:topright)

# Support parameters
support_height = 15.0
support_A = 0.5   # Smaller cross-section for supports
support_I = 0.2   # Smaller moment of inertia for supports
n_elem_per_support = 8

supports = create_uniform_supports(1, L, n_node, support_height, n_elem_per_support, support_A, support_I)

# Bridge with supports
bws = BridgeWithSupports(bo, supports)
M, K = assemble_matrices_with_supports(bws, 20.0)
λs = real(eigvals(Matrix(K), Matrix(M)))  # Convert to dense matrices
ωs = sqrt.(λs[λs .> 1e-6]) ./ (2π)  # Convert to Hz, filter out zero modes

println("Nodes: $(bws.total_nodes), DOFs: $(bws.total_dofs)")
println("Support locations (bridge nodes): $(getfield.(supports, :bridge_node))")

# Plot bridge with supports
p = plot_bridge_with_supports(bws)