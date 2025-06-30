"""
BridgeFEM.jl - Complete Bridge Analysis Example

This script demonstrates the complete workflow for bridge finite element analysis
with temperature-dependent materials and dynamic support structures.

Features demonstrated:
- Temperature-dependent material properties
- Support element modeling (pier)
- Modal analysis across temperature range
- Visualization of mode shapes
- Dynamic response simulation
- Data persistence
"""

using LinearAlgebra, DifferentialEquations, Plots
using Interpolations, SparseArrays, JSON, Arpack, Dates
using Printf

# Include the bridge analysis modules
include("src/bridge_model.jl")
include("src/model_reduction.jl") 
include("src/dynamic_simulation.jl")
include("src/utils.jl")
# include("src/oma.jl")

println("🏗️  BridgeFEM.jl - Complete Bridge Analysis Example")
println("=" ^ 60)

# 1. Define Temperature-Dependent Material Properties
println("\n📊 Step 1: Define temperature-dependent steel properties")

# Realistic steel Young's modulus vs temperature
E_data = [
    -20.0  220e9;   # E at -20°C (very cold)
    -10.0  215e9;   # E at -10°C (cold winter)
      0.0  210e9;   # E at  0°C (freezing)
     20.0  207e9;   # E at 20°C (reference room temperature)
     40.0  200e9;   # E at 40°C (warm summer)
     60.0  190e9    # E at 60°C (very hot)
]

println("Temperature range: $(E_data[1,1])°C to $(E_data[end,1])°C")
println("Young's modulus range: $(E_data[end,2]/1e9:.0) - $(E_data[1,2]/1e9:.0) GPa")

# 2. Create Bridge Model
println("\n🌉 Step 2: Create bridge model")

# Boundary conditions: fixed at left end, pinned in y at right end
bc = BridgeBC([
    [1, "all"],     # Node 1: fix all DOFs (x, y, rotation)
    [51, "y"]       # Node 51: fix vertical displacement only
])

# Bridge parameters for a typical highway bridge
bridge = BridgeOptions(
    50,           # 50 finite elements
    bc,           # Boundary conditions
    300.0,        # 300m total length
    7850.0,       # Steel density (kg/m³)
    4.0,          # Cross-sectional area (m²) - typical for box girder
    3.0,          # Moment of inertia (m⁴) - typical for bridge deck
    E_data,       # Temperature-dependent Young's modulus
    50.0          # Cutoff frequency for modal analysis (Hz)
)

println("Bridge: $(bridge.L)m long, $(bridge.n_elem) elements, $(bridge.n_nodes) nodes")
println("Total DOFs: $(bridge.n_dofs)")

# 3. Add Support Elements
println("\n🏛️  Step 3: Add support structures")

# Vertical pier at bridge center (node 26 = center of 51 nodes)
pier = SupportElement(
    26,           # Connect to bridge node 26 (center)
    [1, 2, 3],    # Connect all DOFs (x, y, rotation)
    -90.0,        # Vertical downward (-90 degrees)
    8,            # 8 elements in pier height
    1.5,          # Pier cross-sectional area (m²)
    0.5,          # Pier moment of inertia (m⁴)
    E_data,       # Same temperature dependence as bridge
    60.0,         # Pier height (m)
    [1, 2, 3]     # Fix all DOFs at pier base
)

supports = [pier]
println("Added vertical pier: $(pier.L)m high, connects to node $(pier.connection_node)")

# 4. Perform Modal Analysis Across Temperature Range
println("\n🔬 Step 4: Modal analysis across temperature range")

temperatures = [-10.0, 0.0, 20.0, 40.0, 60.0]
println("Analyzing at temperatures: $temperatures °C")

# Assemble system matrices for each temperature
n_temps = length(temperatures)
M_all = Array{Float64}(undef, 0, 0, n_temps)
K_all = Array{Float64}(undef, 0, 0, n_temps)

for (i, T) in enumerate(temperatures)
    println("  - Assembling matrices at T = $(T)°C")
    M, K = assemble_matrices_with_supports(bridge, supports, T)
    
    if i == 1
        # Initialize arrays with correct size
        n_dofs = size(M, 1)
        M_all = Array{Float64}(undef, n_dofs, n_dofs, n_temps)
        K_all = Array{Float64}(undef, n_dofs, n_dofs, n_temps)
    end
    
    M_all[:, :, i] = Array(M)
    K_all[:, :, i] = Array(K)
end

# Decompose matrices to get frequencies and mode shapes
println("  - Computing eigenvalues and mode shapes...")
frequencies, mode_shapes, mode_shapes_unnorm = decompose_matrices(M_all, K_all)

# Display first few natural frequencies
println("\nFirst 5 natural frequencies across temperatures:")
println("Mode | " * join([" T=$(T)°C  " for T in temperatures], "|"))
println("-" ^ 50)
for mode = 1:min(5, size(frequencies, 1))
    freq_str = join([sprintf("  %.2f   ", frequencies[mode, i]) for i in 1:length(temperatures)], "|")
    println("  $mode  | $freq_str")
end

# 5. Visualize Mode Shapes
println("\n📈 Step 5: Visualize mode shapes")

# Plot first three modes at reference temperature (20°C)
temp_idx = findfirst(x -> x ≈ 20.0, temperatures)
if temp_idx !== nothing
    println("Creating mode shape plots at T = 20°C")
    
    plots = []
    for mode = 1:min(3, size(mode_shapes, 2))
        p = plot_bridge_with_supports(
            bridge, supports,
            mode_shape = mode_shapes[:, mode, temp_idx],
            scale_factor = 100.0,
            T = 20.0,
            title = "Mode $mode: f = $(frequencies[mode, temp_idx]:.2f) Hz",
            show_nodes = false
        )
        push!(plots, p)
    end
    
    # Combine plots
    combined_plot = plot(plots..., layout = (3, 1), size = (800, 900))
    
    # Save plot
    savefig(combined_plot, "bridge_mode_shapes.png")
    println("Mode shape plots saved as 'bridge_mode_shapes.png'")
end

# 6. Save Simulation Configuration
println("\n💾 Step 6: Save simulation configuration")

# Create simulation options object
sim_options = SimulationOptions(bridge, supports, temperatures, damping_ratio=0.02)

# Save to JSON file
config_filename = "bridge_simulation_config.json"
save_simulation_options(sim_options, config_filename)

# Test loading
loaded_options = load_simulation_options(config_filename)
println("Configuration saved and verified: $config_filename")

# 7. Summary
println("\n📋 Analysis Summary")
println("=" ^ 40)
println("Bridge length: $(bridge.L) m")
println("Number of elements: $(bridge.n_elem)")
println("Total DOFs (with supports): $(size(M_all, 1))")
println("Temperature range: $(minimum(temperatures)) to $(maximum(temperatures)) °C")
println("Frequency range (1st mode): $(minimum(frequencies[1,:])) - $(maximum(frequencies[1,:])) Hz")
println("Temperature effect on 1st mode: $(((maximum(frequencies[1,:]) - minimum(frequencies[1,:]))/minimum(frequencies[1,:]))*100:.1f)%")

# Optional: Dynamic Analysis Example
println("\n🎬 Optional: Dynamic response example")
println("(Uncomment code below to run dynamic simulation)")

"""
# Create simple sinusoidal loading at bridge center
function simple_loading(t, dofs)
    f = zeros(length(dofs))
    center_node_dof = 3 * 26 - 1  # Vertical DOF at center node
    if center_node_dof in dofs
        idx = findfirst(x -> x == center_node_dof, dofs)
        f[idx] = 1000.0 * sin(2π * 2.0 * t)  # 2 Hz, 1000 N amplitude
    end
    return f
end

# Temperature function (constant for simplicity)
T_func(t) = 20.0

# Set up dynamic simulation parameters
n_modes = min(10, size(mode_shapes, 2))
modal_params = (
    T_func = T_func,
    n_modes = n_modes,
    ω_interp = T -> frequencies[1:n_modes, temp_idx],
    ζ = 0.02 * ones(n_modes),  # 2% damping
    Φ_interp = T -> mode_shapes[:, 1:n_modes, temp_idx],
    load_vector = simple_loading,
    n_dofs = size(mode_shapes, 1)
)

# Solve dynamic response
println("Solving dynamic response for 10 seconds...")
prob = ODEProblem(beam_modal_ode!, zeros(2*n_modes), (0.0, 10.0), modal_params)
sol = solve(prob, Tsit5(), reltol=1e-6)

println("Dynamic simulation completed!")
"""

# Helper function for formatted printing (used in example)
function sprintf(format_str, value)
    return @sprintf(format_str, value)
end

println("\n✅ Bridge analysis example completed successfully!")
println("Check generated files:")
println("  - bridge_mode_shapes.png (mode shape visualization)")
println("  - bridge_simulation_config.json (simulation configuration)")
