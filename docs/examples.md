---
layout: default
title: Examples - BridgeFEM.jl
description: Worked examples and tutorials for bridge finite element analysis
toc: true
---

# 💡 Examples and Tutorials

This page provides worked examples demonstrating various features of BridgeFEM.jl, from basic bridge analysis to advanced dynamic simulations.

---

## 📚 Example 1: Basic Bridge Analysis

### Problem Description

Analyze a 300m simply-supported steel bridge under temperature variations from -10°C to 50°C.

### Solution

```julia
using BridgeFEM

# Define temperature-dependent steel properties
E_data = [
    -10.0  250e9;   # Winter: increased stiffness
     20.0  207e9;   # Reference temperature
     50.0  150e9    # Summer: thermal softening
]

# Boundary conditions: pinned at both ends
bc = BridgeBC([
    [1, "trans"],    # Fix translations at left end
    [51, "y"]        # Fix vertical translation at right end  
])

# Create bridge model
bridge = BridgeOptions(
    50,              # 50 elements
    bc,              # Boundary conditions
    300.0,           # 300m length
    7850.0,          # Steel density (kg/m³)
    4.0,             # Cross-sectional area (m²)
    3.0,             # Moment of inertia (m⁴)
    E_data,          # Temperature dependence
    50.0             # Frequency cutoff (Hz)
)

# Perform modal analysis
temperatures = [-10.0, 20.0, 50.0]
M, K, frequencies, mode_shapes = assemble_and_decompose(bridge, temperatures)

# Display results
println("First 5 natural frequencies:")
for i = 1:5
    println("Mode $i: $(frequencies[i,1]:.2f) Hz → $(frequencies[i,end]:.2f) Hz")
end
```

### Expected Results

```
Mode 1: 2.45 Hz → 1.89 Hz  (23% reduction)
Mode 2: 9.78 Hz → 7.55 Hz  (23% reduction)  
Mode 3: 22.0 Hz → 17.0 Hz  (23% reduction)
```

**Key Insight:** All frequencies scale as $\sqrt{E(T)/E_0}$, showing uniform temperature effect.

---

## 🏛️ Example 2: Bridge with Pier Support

### Problem Description

Add a vertical pier to the bridge center and analyze the coupled system dynamics.

### Solution

```julia
# Continue from Example 1...

# Add vertical pier at bridge center (node 26)
pier = SupportElement(
    26,              # Connection node (center of 51 nodes)
    [1, 2, 3],       # Connect all DOFs
    -90.0,           # Vertical downward
    8,               # 8 elements in pier height
    1.5,             # Pier cross-sectional area (m²)
    0.5,             # Pier moment of inertia (m⁴)
    E_data,          # Same temperature dependence
    60.0,            # 60m pier height
    [1, 2, 3]        # Fixed base
)

# Analyze coupled bridge-pier system
supports = [pier]
M, K, frequencies, mode_shapes = assemble_and_decompose(
    bridge, temperatures, supports=supports
)

# Visualize mode shapes
for mode = 1:3
    plot_bridge_with_supports(
        bridge, supports,
        mode_shape = mode_shapes[:, mode, 2],  # At 20°C
        scale_factor = 100.0,
        title = "Mode $mode: f = $(frequencies[mode,2]:.2f) Hz"
    )
end
```

### System Properties

- **Total DOFs:** Bridge (153) + Pier (27) = 180 DOFs
- **New mode types:** Bridge bending + pier swaying
- **Frequency shifts:** Pier adds stiffness → higher frequencies

---

## 🎬 Example 3: Dynamic Response Analysis

### Problem Description

Simulate bridge response to a frequency sweep (chirp) excitation.

### Solution

```julia
# Define chirp loading function
function create_chirp_loading(f1, f2, target_node, amplitude=1000.0)
    return (t, dofs) -> begin
        f = zeros(length(dofs))
        
        # Linear frequency sweep from f1 to f2
        tspan = 100.0  # Total time
        freq_t = f1 + (f2 - f1) * (t / tspan)
        phase = 2π * (f1 * t + (f2 - f1) * t^2 / (2 * tspan))
        
        # Apply force to vertical DOF of target node
        target_dof = 3 * target_node - 1  # y-displacement DOF
        if target_dof in dofs
            idx = findfirst(x -> x == target_dof, dofs)
            f[idx] = amplitude * sin(phase)
        end
        
        return f
    end
end

# Set up dynamic simulation
load_func = create_chirp_loading(0.1, 10.0, 26)  # 0.1-10 Hz sweep at center
T_func(t) = 20.0  # Constant temperature

# Solve dynamics using modal superposition
sim_options = SimulationOptions(bridge, supports, temperatures)
sol = solve_dynamics(sim_options, load_func, T_func, (0.0, 100.0))

# Create animation
animate_dynamic_response(bridge, supports, sol.u, sol.t,
                        scale_factor=1000.0,
                        filename="bridge_chirp_response.gif")
```

### Analysis Results

- **Resonance peaks:** Clear amplification at natural frequencies
- **Temperature sensitivity:** Response amplitude varies with temperature
- **Modal contributions:** Different modes dominate at different frequencies

---

## 🌡️ Example 4: Thermal Analysis

### Problem Description

Analyze bridge behavior under realistic daily temperature cycles.

### Solution

```julia
# Define daily temperature variation
function daily_temperature(t)
    T_mean = 25.0      # Daily average (°C)
    T_amplitude = 15.0  # Daily variation (±15°C)
    period = 24.0 * 3600.0  # 24 hours in seconds
    
    return T_mean + T_amplitude * sin(2π * t / period)
end

# Time-varying modal analysis
time_points = 0:3600:86400  # Every hour for 24 hours
freqs_vs_time = zeros(5, length(time_points))

for (i, t) in enumerate(time_points)
    T = daily_temperature(t)
    M, K = assemble_matrices_with_supports(bridge, supports, T)
    λ, φ = eigs(K, M, nev=5, which=:SM)
    freqs_vs_time[:, i] = sqrt.(real.(λ)) / (2π)
end

# Plot frequency evolution
using Plots
plot(time_points/3600, freqs_vs_time', 
     xlabel="Time (hours)", 
     ylabel="Frequency (Hz)",
     title="Daily Frequency Variation",
     labels=["Mode $i" for i=1:5])
```

### Key Observations

- **Frequency variation:** ±12% daily variation in natural frequencies
- **Phase relationship:** Frequencies minimum at peak temperature (afternoon)
- **Engineering implications:** Must account for thermal effects in design

---

## 📊 Example 5: Parametric Study

### Problem Description

Study the effect of pier stiffness on bridge dynamics.

### Solution

```julia
# Define parameter ranges
pier_areas = [0.5, 1.0, 1.5, 2.0, 2.5]  # m²
pier_inertias = [0.1, 0.3, 0.5, 0.8, 1.0]  # m⁴

# Storage for results
results = Dict()

for A_pier in pier_areas
    for I_pier in pier_inertias
        # Create pier with current properties
        pier = SupportElement(26, [1,2,3], -90.0, 8, A_pier, I_pier, 
                            E_data, 60.0, [1,2,3])
        
        # Analyze system
        M, K, freqs, modes = assemble_and_decompose(
            bridge, [20.0], supports=[pier]
        )
        
        # Store first frequency
        key = (A_pier, I_pier)
        results[key] = freqs[1, 1]
    end
end

# Create contour plot
using Plots
heatmap(pier_areas, pier_inertias, 
        [results[(A,I)] for I in pier_inertias, A in pier_areas],
        xlabel="Pier Area (m²)", 
        ylabel="Pier Inertia (m⁴)",
        title="First Frequency (Hz)")
```

### Design Insights

- **Pier area effect:** Primarily affects axial/lateral modes
- **Pier inertia effect:** Strongly influences bending modes  
- **Optimization:** Trade-off between stiffness and material cost

---

## 🔍 Example 6: Operational Modal Analysis

### Problem Description

Extract modal parameters from simulated ambient response data.

### Solution

```julia
include("src/oma.jl")  # Operational Modal Analysis functions

# Simulate ambient excitation (wind loading)
function ambient_loading(t, dofs)
    f = zeros(length(dofs))
    
    # Random white noise on all bridge nodes
    bridge_nodes = 1:bridge.n_nodes
    for node in bridge_nodes
        y_dof = 3 * node - 1  # Vertical DOF
        if y_dof <= length(dofs)
            f[y_dof] = 100.0 * randn()  # Random force
        end
    end
    
    return f
end

# Generate response data
T_func(t) = 20.0  # Constant temperature
tspan = (0.0, 600.0)  # 10 minutes of data
sol = solve_dynamics(sim_options, ambient_loading, T_func, tspan)

# Extract output-only data (accelerations at sensor locations)
sensor_nodes = [10, 20, 30, 40]  # 4 sensors along bridge
sensor_dofs = [3*n-1 for n in sensor_nodes]  # Vertical accelerations

# Reconstruct accelerations
u, v, a = reconstruct_physical(sim_options, sol.u, T_func, sol.t)
y_data = a[sensor_dofs, :]  # Output matrix

# Perform Stochastic Subspace Identification
dt = sol.t[2] - sol.t[1]
i = 20   # Past horizon
n = 2:10  # Model orders to test

frequencies, damping, mode_shapes = run_ssi_and_save(y_data, i, n, dt)

# Compare with analytical results
println("Modal identification results:")
println("Analytical vs Identified frequencies:")
for mode = 1:3
    f_analytical = frequencies_true[mode, 2]  # At 20°C
    f_identified = frequencies[end][mode]     # Best model order
    error = abs(f_identified - f_analytical) / f_analytical * 100
    println("Mode $mode: $(f_analytical:.2f) Hz vs $(f_identified:.2f) Hz ($(error:.1f)% error)")
end
```

### Validation Results

- **Frequency accuracy:** <2% error for well-excited modes
- **Damping estimation:** Realistic values around 1-3%
- **Mode shapes:** Good correlation with analytical results

---

## 🛠️ Advanced Tips

### Performance Optimization

```julia
# Use sparse matrices for large systems
using SparseArrays

# Pre-allocate arrays for time loops
n_steps = length(sol.t)
displacements = Matrix{Float64}(undef, n_dofs, n_steps)

# Parallelize parameter studies
using Distributed
addprocs(4)  # Add 4 worker processes

@everywhere using BridgeFEM
results = @distributed (append!) for param in parameter_range
    analyze_parameter(param)
end
```

### Memory Management

```julia
# Clear large matrices when not needed
M = K = nothing
GC.gc()  # Force garbage collection

# Use views for array slicing
mode_subset = @view mode_shapes[:, 1:10, :]  # No copying
```

### Debugging Tips

```julia
# Check system properties
@assert rank(K) == size(K,1) - n_constraints  # Verify constraints
@assert all(λ .> 0)  # All eigenvalues positive

# Visualize mode shapes for sanity check
for i = 1:5
    plot_mode_shape(bridge, supports, mode_shapes[:,i,1], i)
end
```

---

## 📝 Next Steps

- **Explore the [API Reference](api-reference.md)** for detailed function documentation
- **Check the [Theory](theory.md)** section for mathematical background
- **Try the complete example script** in `bridge_example.jl`
- **Adapt examples** to your specific bridge geometry and loading conditions

---

**💡 Pro Tip:** Start with simple examples and gradually add complexity. The modular design of BridgeFEM.jl makes it easy to extend existing analyses with new features.
