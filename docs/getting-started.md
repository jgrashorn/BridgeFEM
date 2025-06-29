---
layout: default
title: Getting Started - BridgeFEM.jl
---

# Getting Started

This guide will help you install BridgeFEM.jl and run your first simulation.

## Installation

### Requirements

- Julia 1.6 or higher
- Recommended: 8GB RAM, multi-core processor for large models

### Install BridgeFEM.jl

```julia
using Pkg

# Option 1: Install from GitHub (development version)
Pkg.add(url="https://github.com/jgrashorn/BridgeFEM.jl")

# Option 2: Development installation
Pkg.develop(url="https://github.com/jgrashorn/BridgeFEM.jl")
```

### Install Dependencies

```julia
# Required packages (installed automatically)
using LinearAlgebra, DifferentialEquations, Plots
using Interpolations, SparseArrays, JSON, Arpack, Dates
```

### Test Installation

```julia
using BridgeFEM

# Create a simple bridge
bc = BridgeBC([[1, "all"], [11, "y"]])  
bridge = BridgeOptions(10, bc, 100.0, 7800.0, 1.0, 0.1, [20.0 207e9], 20.0)

# Run modal analysis
M, K, frequencies, modes = assemble_and_decompose(bridge, [20.0])

println("First 3 natural frequencies: ", frequencies[1:3, 1], " Hz")
```

If this runs without error, you're ready to go!

---

## Basic Concepts

### Bridge Definition

A bridge is defined by the `BridgeOptions` struct:

```julia
# Material properties with temperature dependence
E_data = [
    -10.0  250e9;   # Young's modulus at -10°C
     20.0  207e9;   # Young's modulus at 20°C
     50.0  150e9    # Young's modulus at 50°C
]

# Boundary conditions: fix all DOFs at node 1, y-displacement at last node
bc = BridgeBC([
    [1, "all"],      # x, y, rotation fixed
    [21, "y"]        # only y-displacement fixed
])

# Create bridge: 20 elements, 200m long, steel properties
bridge = BridgeOptions(
    20,              # Number of elements
    bc,              # Boundary conditions
    200.0,           # Length (m)
    7800.0,          # Density (kg/m³)
    2.0,             # Cross-sectional area (m²)
    0.5,             # Moment of inertia (m⁴)
    E_data,          # Temperature-dependent Young's modulus
    50.0             # Cutoff frequency for modes (Hz)
)
```

### Support Elements

Add piers, cables, or other supports:

```julia
# Vertical pier at bridge center
pier = SupportElement(
    11,              # Connect to node 11 (middle of bridge)
    [1, 2, 3],       # Connect x, y, rotation DOFs
    -90.0,           # Angle: -90° = vertical downward
    5,               # 5 elements in pier
    1.0,             # Cross-sectional area (m²)
    0.1,             # Moment of inertia (m⁴)
    E_data,          # Same material as bridge
    30.0,            # Pier height (m)
    [1, 2, 3]        # Fix all DOFs at pier base
)
```

---

## Your First Simulation

### 1. Modal Analysis

```julia
using BridgeFEM, Plots

# Define bridge and pier (as above)
# ...

# Perform modal analysis at different temperatures
temperatures = [-10.0, 0.0, 20.0, 40.0, 50.0]
M, K, frequencies, mode_shapes = assemble_and_decompose(bridge, temperatures, supports=[pier])

# Plot frequency vs temperature for first 3 modes
plot(temperatures, frequencies[1:3, :]', 
     xlabel="Temperature (°C)", 
     ylabel="Frequency (Hz)",
     labels=["Mode 1" "Mode 2" "Mode 3"],
     linewidth=2, marker=:circle)
title!("Natural Frequencies vs Temperature")
```

### 2. Visualize Mode Shapes

```julia
# Plot first mode shape at 20°C
mode1_20C = mode_shapes[:, 1, 3]  # Mode 1, temperature index 3 (20°C)

plot_bridge_with_supports(bridge, [pier], 
                         mode_shape=mode1_20C, 
                         scale_factor=100.0,
                         title="First Mode Shape at 20°C")
```

### 3. Dynamic Response

```julia
# Create frequency sweep loading function
function create_chirp_loading(f1, f2, magnitude, target_node)
    return (t, dof) -> begin
        f = zeros(length(dof))
        
        # Linear frequency sweep
        freq_t = f1 + (f2 - f1) * (t / 10.0)  # 10 second sweep
        phase = 2π * (f1 * t + (f2 - f1) * t^2 / 20.0)
        
        # Apply to y-DOF of target node
        y_dof = 3 * (target_node - 1) + 2
        f[y_dof] = magnitude * sin(phase)
        
        return f
    end
end

# Set up dynamic simulation
load_vector = create_chirp_loading(0.1, 10.0, 10000.0, 11)  # Sweep 0.1-10 Hz, 10kN force

# Create interpolators for temperature-dependent properties
using Interpolations
ω_interp = interpolate((1:size(frequencies,1), temperatures), frequencies, Gridded(Linear()))
Φ_interp = interpolate((1:size(mode_shapes,1), 1:size(mode_shapes,2), temperatures), mode_shapes, Gridded(Linear()))

# Temperature function (constant for this example)
T_func = t -> 20.0

# Solve dynamic system
n_modes = size(frequencies, 1)
u0 = zeros(2 * n_modes)  # Initial conditions: [displacements; velocities]

params = (
    n_modes = n_modes,
    n_dofs = size(mode_shapes, 1),
    T_func = T_func,
    ω_interp = (modes, T) -> 2π .* ω_interp(modes, T),  # Convert Hz to rad/s
    ζ = fill(0.02, n_modes),  # 2% damping for all modes
    Φ_interp = Φ_interp,
    load_vector = load_vector
)

# Solve ODE
using DifferentialEquations
prob = ODEProblem(beam_modal_ode!, u0, (0.0, 10.0), params)
sol = solve(prob, BS3(), saveat=0.02)

# Reconstruct physical displacements
q = reduce(hcat, sol.u)
u_phys, v_phys = reconstruct_physical(bridge, q, Φ_interp, T_func, sol.t, supports=[pier])

# Plot response at bridge center
center_node = 11
y_dof = 3 * (center_node - 1) + 2
plot(sol.t, u_phys[y_dof, :], 
     xlabel="Time (s)", ylabel="Displacement (m)",
     title="Bridge Center Response to Frequency Sweep",
     linewidth=2)
```

### 4. Create Animation

```julia
# Animate the dynamic response
anim = animate_dynamic_response(bridge, [pier], u_phys, sol.t,
                              scale_factor=10000.0,  # Amplify for visibility
                              n_frames=100,
                              fps=20,
                              filename="bridge_frequency_sweep.gif")

# The animation is automatically saved as a GIF
```

---

## Configuration Management

### Save Simulation Setup

```julia
# Create complete simulation configuration
sim_opts = SimulationOptions(bridge, [pier], temperatures, damping_ratio=0.02)

# Save to JSON for reproducibility
save_simulation_options(sim_opts, "my_bridge_config.json")
```

### Load and Reuse Configuration

```julia
# Load configuration from file
loaded_opts = load_simulation_options("my_bridge_config.json")

# Extract components
bridge_reloaded = loaded_opts.bridge
supports_reloaded = loaded_opts.supports
temps_reloaded = loaded_opts.temperatures

# Run analysis with loaded configuration
M, K, frequencies, mode_shapes = assemble_and_decompose(
    bridge_reloaded, temps_reloaded, supports=supports_reloaded)
```

---

## Common Patterns

### Parameter Studies

```julia
# Study effect of pier stiffness
pier_areas = [0.5, 1.0, 2.0, 4.0]
first_frequencies = zeros(length(pier_areas))

for (i, A) in enumerate(pier_areas)
    # Create pier with different stiffness
    pier_var = SupportElement(11, [1,2,3], -90.0, 5, A, A^2/12, E_data, 30.0, [1,2,3])
    
    # Analyze
    _, _, freqs, _ = assemble_and_decompose(bridge, [20.0], supports=[pier_var])
    first_frequencies[i] = freqs[1, 1]
end

# Plot results
plot(pier_areas, first_frequencies,
     xlabel="Pier Cross-sectional Area (m²)",
     ylabel="First Natural Frequency (Hz)",
     marker=:circle, linewidth=2)
```

### Multiple Supports

```julia
# Bridge with multiple piers
pier1 = SupportElement(8, [1,2,3], -90.0, 4, 1.0, 0.1, E_data, 25.0, [1,2,3])
pier2 = SupportElement(11, [1,2,3], -90.0, 5, 1.2, 0.12, E_data, 30.0, [1,2,3])  
pier3 = SupportElement(14, [1,2,3], -90.0, 4, 1.0, 0.1, E_data, 25.0, [1,2,3])

supports = [pier1, pier2, pier3]

# Analysis with multiple supports
M, K, frequencies, mode_shapes = assemble_and_decompose(bridge, temperatures, supports=supports)

# Visualization
plot_bridge_with_supports(bridge, supports, 
                         mode_shape=mode_shapes[:, 1, 3],
                         scale_factor=50.0,
                         title="Bridge with Three Piers")
```

### Different Boundary Conditions

```julia
# Simply supported bridge
bc_simply = BridgeBC([[1, "trans"], [21, "y"]])  # Pin and roller

# Cantilever bridge  
bc_cantilever = BridgeBC([[1, "all"]])  # Fixed at one end

# Continuous bridge (no support reactions, only pier constraints)
bc_continuous = BridgeBC([])  # Free ends, supported by piers only
```

---

## Performance Tips

### Large Models
For models with many DOFs (>1000):

```julia
# Use higher cutoff frequency to capture more modes
bridge.cutoff_freq = 100.0

# Reduce mode count for faster dynamics
n_modes_keep = 20
frequencies_reduced = frequencies[1:n_modes_keep, :]
mode_shapes_reduced = mode_shapes[:, 1:n_modes_keep, :]
```

### Long Simulations
For time simulations longer than 100 seconds:

```julia
# Use larger time steps
prob = ODEProblem(beam_modal_ode!, u0, (0.0, 1000.0), params)
sol = solve(prob, BS3(), saveat=0.1)  # 0.1s instead of 0.02s

# Or use manual integration for better control
include("manual_integrators.jl")
t_manual, u_manual = manual_rk4_integration(u0, (0.0, 1000.0), 0.05, params)
```

### Memory Management
For animations and large datasets:

```julia
# Subsample for animation
time_subsample = sol.t[1:5:end]  # Every 5th time point
u_subsample = u_phys[:, 1:5:end]

animate_dynamic_response(bridge, supports, u_subsample, time_subsample,
                        n_frames=50)  # Fewer frames
```

---

## Troubleshooting

### Common Issues

**1. Negative eigenvalues**
```julia
# Check boundary conditions
println("Number of constrained DOFs: ", length(bc.conds))

# Verify support connectivity
supports_dof_maps, total_dofs = create_support_dof_mapping(bridge, supports)
println("Total DOFs: ", total_dofs)
```

**2. Integration instability**
```julia
# Reduce time step
sol = solve(prob, BS3(), saveat=0.001)  # Smaller time step

# Or use implicit method
sol = solve(prob, Rodas4(), saveat=0.01)  # Implicit solver
```

**3. Memory issues**
```julia
# Monitor memory usage
using Profile
@profile solve(prob, BS3(), saveat=0.01)
```

### Getting Help

- **Documentation**: See [API Reference](api-reference.html) for detailed function descriptions
- **Examples**: Check [Examples](examples.html) for more complex scenarios
- **Issues**: Report bugs on [GitHub Issues](https://github.com/jgrashorn/BridgeFEM.jl/issues)

---

**Next**: [Theory](theory.html) | **See also**: [Examples](examples.html)
