---
layout: default
title: API Reference - BridgeFEM.jl
---

# API Reference

Complete reference for all functions and types in BridgeFEM.jl.

## Core Types

### BridgeOptions

```julia
BridgeOptions(n_elem, bc_nodes, L, ρ, A, I, E_T, cutoff_freq)
```

Main structure defining a bridge finite element model.

**Arguments:**
- `n_elem::Int`: Number of finite elements
- `bc_nodes::BridgeBC`: Boundary conditions  
- `L::Float64`: Bridge length (m)
- `ρ::Float64`: Material density (kg/m³)
- `A::Float64`: Cross-sectional area (m²)
- `I::Float64`: Moment of inertia (m⁴)
- `E_T::Matrix{Float64}`: Temperature-Young's modulus data `[T E; ...]`
- `cutoff_freq::Float64`: Maximum frequency for modal analysis (Hz)

**Fields:**
- `n_nodes::Int`: Number of nodes (computed as `n_elem + 1`)
- `n_dofs::Int`: Number of degrees of freedom (3 per node)
- `E::Function`: Young's modulus interpolation function `E(T)`

### SupportElement

```julia
SupportElement(connection_node, connection_dofs, angle, n_elem, A, I, E_T, L, bc_bottom)
```

Auxiliary support structure connected to the main bridge.

**Arguments:**
- `connection_node::Int`: Bridge node to connect to
- `connection_dofs::Vector{Int}`: DOFs to connect `[1,2,3] = [x,y,θ]`
- `angle::Float64`: Orientation angle in degrees
- `n_elem::Int`: Number of elements in support
- `A::Float64`: Cross-sectional area (m²)
- `I::Float64`: Moment of inertia (m⁴)
- `E_T::Matrix{Float64}`: Temperature-dependent Young's modulus data
- `L::Float64`: Support length (m)
- `bc_bottom::Vector{Int}`: Constrained DOFs at support base

**Coordinate Convention:**
- 0° = horizontal right, -90° = vertical down
- First node connects to bridge, last node is constrained

### BridgeBC

```julia
BridgeBC(conditions)
```

Boundary condition specification for bridge nodes.

**Arguments:**
- `conditions::Vector{Vector{Any}}`: List of `[node_number, constraint_type]`

**Constraint Types:**
- String: `"all"`, `"trans"`, `"x"`, `"y"`, `"ϕ"`
- Vector: `[1, 2, 3]` for explicit DOF specification

### SimulationOptions

```julia
SimulationOptions(bridge, supports, temperatures; damping_ratio=0.02)
```

Complete simulation configuration for data persistence.

**Arguments:**
- `bridge::BridgeOptions`: Bridge configuration
- `supports::Vector{SupportElement}`: Support structures
- `temperatures::Vector{Float64}`: Temperature range
- `damping_ratio::Float64`: Modal damping ratio

---

## Assembly Functions

### assemble_matrices

```julia
assemble_matrices(bridge, T=20.0) -> (M, K)
```

Assemble global mass and stiffness matrices for bridge only.

**Arguments:**
- `bridge::BridgeOptions`: Bridge configuration
- `T::Float64`: Temperature (°C)

**Returns:**
- `M::SparseMatrixCSC`: Global mass matrix
- `K::SparseMatrixCSC`: Global stiffness matrix

### assemble_matrices_with_supports

```julia
assemble_matrices_with_supports(bridge, supports, T=20.0) -> (M, K)
```

Assemble expanded system matrices including supports.

**Arguments:**
- `bridge::BridgeOptions`: Bridge configuration
- `supports::Vector{SupportElement}`: Support structures  
- `T::Float64`: Temperature (°C)

**Returns:**
- `M::SparseMatrixCSC`: Expanded mass matrix
- `K::SparseMatrixCSC`: Expanded stiffness matrix

### assemble_and_decompose

```julia
assemble_and_decompose(bridge, temperatures; supports=[]) -> (M, K, λs, vectors, vectors_unnorm)
```

Complete workflow: assembly + eigenvalue decomposition.

**Arguments:**
- `bridge::BridgeOptions`: Bridge configuration
- `temperatures::Vector{Float64}`: Temperature points
- `supports::Vector{SupportElement}`: Support structures (optional)

**Returns:**
- `M::Array{Float64,3}`: Mass matrices `[dof, dof, temperature]`
- `K::Array{Float64,3}`: Stiffness matrices `[dof, dof, temperature]`
- `λs::Matrix{Float64}`: Natural frequencies `[mode, temperature]` (Hz)
- `vectors::Array{Float64,3}`: Normalized mode shapes `[dof, mode, temperature]`
- `vectors_unnorm::Array{Float64,3}`: Unnormalized mode shapes

---

## Analysis Functions

### decompose_matrices

```julia
decompose_matrices(M, K) -> (eigenvalues, eigenvectors)
```

Solve generalized eigenvalue problem for modal analysis.

**Arguments:**
- `M::AbstractMatrix`: Mass matrix
- `K::AbstractMatrix`: Stiffness matrix

**Returns:**
- `eigenvalues::Vector{Float64}`: Natural frequencies (rad²/s²)
- `eigenvectors::Matrix{Float64}`: Mode shapes (normalized)

### beam_modal_ode!

```julia
beam_modal_ode!(du, u, params, t)
```

ODE function for modal coordinate integration.

**Arguments:**
- `du::Vector{Float64}`: Output derivative vector
- `u::Vector{Float64}`: Current state `[modal_disp; modal_vel]`
- `params::NamedTuple`: System parameters
- `t::Float64`: Current time

**Parameters Structure:**
```julia
params = (
    n_modes = n_modes,
    n_dofs = total_dofs,
    T_func = t -> T(t),           # Temperature function
    ω_interp = ω_interpolator,    # Frequency interpolator
    ζ = damping_ratios,           # Modal damping ratios
    Φ_interp = mode_interpolator, # Mode shape interpolator
    load_vector = force_function  # External loading function
)
```

### reconstruct_physical

```julia
reconstruct_physical(bridge, q_full, Φ_interp, T_func, time; supports=[]) -> (u, du)
```

Convert modal coordinates to physical displacements.

**Arguments:**
- `bridge::BridgeOptions`: Bridge configuration
- `q_full::Matrix{Float64}`: Modal coordinate history `[modal_dof, time]`
- `Φ_interp`: Mode shape interpolator
- `T_func::Function`: Temperature function `T(t)`
- `time::Vector{Float64}`: Time vector
- `supports::Vector{SupportElement}`: Support structures (optional)

**Returns:**
- `u::Matrix{Float64}`: Physical displacements `[dof, time]`
- `du::Matrix{Float64}`: Physical velocities `[dof, time]`

---

## Element Functions

### frame_elem_stiffness

```julia
frame_elem_stiffness(EA, EI, L_e) -> k
```

Compute 6×6 stiffness matrix for 2D frame element.

**Arguments:**
- `EA::Float64`: Axial stiffness (N)
- `EI::Float64`: Bending stiffness (N⋅m²)
- `L_e::Float64`: Element length (m)

**Returns:**
- `k::Matrix{Float64}`: 6×6 element stiffness matrix

**DOF Ordering:** `[u₁, v₁, θ₁, u₂, v₂, θ₂]`

### frame_elem_mass

```julia
frame_elem_mass(ρ, A, L) -> m
```

Compute 6×6 consistent mass matrix for 2D frame element.

**Arguments:**
- `ρ::Float64`: Material density (kg/m³)
- `A::Float64`: Cross-sectional area (m²)
- `L::Float64`: Element length (m)

**Returns:**
- `m::Matrix{Float64}`: 6×6 element mass matrix

---

## Support Functions

### create_support_dof_mapping

```julia
create_support_dof_mapping(bridge, supports) -> (dof_maps, total_dofs)
```

Create DOF mapping between bridge and support systems.

**Arguments:**
- `bridge::BridgeOptions`: Bridge configuration
- `supports::Vector{SupportElement}`: Support structures

**Returns:**
- `dof_maps::Vector{Vector{Int}}`: Mapping from local to global DOFs for each support
- `total_dofs::Int`: Total number of DOFs in expanded system

### create_expanded_transformation

```julia
create_expanded_transformation(angle, n_nodes) -> T_expanded
```

Create transformation matrix for rotating support to global coordinates.

**Arguments:**
- `angle::Float64`: Rotation angle (degrees)
- `n_nodes::Int`: Number of nodes in support

**Returns:**
- `T_expanded::Matrix{Float64}`: Block-diagonal transformation matrix

---

## Visualization Functions

### plot_bridge_with_supports

```julia
plot_bridge_with_supports(bridge, supports; kwargs...) -> plot
```

Plot bridge structure with optional mode shape visualization.

**Arguments:**
- `bridge::BridgeOptions`: Bridge configuration
- `supports::Vector{SupportElement}`: Support structures

**Keyword Arguments:**
- `mode_shape::Vector{Float64}`: Mode shape displacements (optional)
- `scale_factor::Float64=1.0`: Scaling factor for deformation
- `title::String="Bridge Structure"`: Plot title
- `show_nodes::Bool=true`: Whether to show node markers
- `T::Float64=20.0`: Temperature for material properties

**Returns:**
- Plot object for display or further modification

### animate_dynamic_response

```julia
animate_dynamic_response(bridge, supports, u, time; kwargs...) -> animation
```

Create animation of dynamic bridge response.

**Arguments:**
- `bridge::BridgeOptions`: Bridge configuration
- `supports::Vector{SupportElement}`: Support structures
- `u::Matrix{Float64}`: Displacement history `[dof, time]`
- `time::Vector{Float64}`: Time vector

**Keyword Arguments:**
- `scale_factor::Float64=1.0`: Amplification for visualization
- `n_frames::Int=100`: Number of animation frames
- `fps::Int=20`: Frames per second
- `filename::String="bridge_dynamics.gif"`: Output GIF filename
- `show_nodes::Bool=false`: Whether to show node markers
- `title_prefix::String="Bridge Dynamic Response"`: Title prefix

**Returns:**
- Animation object (automatically saved as GIF)

---

## Data Management

### save_simulation_options

```julia
save_simulation_options(opts, filename)
```

Save complete simulation setup to JSON file.

**Arguments:**
- `opts::SimulationOptions`: Simulation configuration
- `filename::String`: Output JSON filename

### load_simulation_options

```julia
load_simulation_options(filename) -> SimulationOptions
```

Load simulation setup from JSON file.

**Arguments:**
- `filename::String`: Input JSON filename

**Returns:**
- `SimulationOptions`: Loaded simulation configuration

### bridge_options_to_dict / dict_to_bridge_options

```julia
bridge_options_to_dict(bridge) -> Dict
dict_to_bridge_options(dict) -> BridgeOptions
```

Convert between `BridgeOptions` and dictionary for JSON serialization.

### support_element_to_dict / dict_to_support_element

```julia
support_element_to_dict(support) -> Dict  
dict_to_support_element(dict) -> SupportElement
```

Convert between `SupportElement` and dictionary for JSON serialization.

---

## Utility Functions

### BCTypes

```julia
BCTypes::Dict{String, Any}
```

Dictionary mapping boundary condition strings to DOF indices:
- `"all"` → `[1,2,3]` (all DOFs)
- `"trans"` → `[1,2]` (translations only)
- `"x"` → `1` (x-translation)
- `"y"` → `2` (y-translation)  
- `"ϕ"` → `3` (rotation)

### transformation_matrix

```julia
transformation_matrix(θ) -> T
```

Create 6×6 transformation matrix for 2D rotation.

**Arguments:**
- `θ::Float64`: Rotation angle (degrees)

**Returns:**
- `T::Matrix{Float64}`: 6×6 transformation matrix

### apply_bc!

```julia
apply_bc!(M, K, bc_dofs)
```

Apply boundary conditions to mass and stiffness matrices (in-place).

**Arguments:**
- `M::AbstractMatrix`: Mass matrix (modified)
- `K::AbstractMatrix`: Stiffness matrix (modified)
- `bc_dofs::Vector`: Boundary condition specification

---

## Loading Functions

Common loading patterns for dynamic analysis:

### Frequency Sweep (Chirp)

```julia
# Linear frequency sweep
load_vector = (t, dof) -> begin
    f = zeros(length(dof))
    f1, f2 = 0.1, 10.0  # Hz
    freq_t = f1 + (f2 - f1) * (t / T_total)
    phase = 2π * (f1 * t + (f2 - f1) * t^2 / (2 * T_total))
    f[target_dof] = magnitude * sin(phase)
    return f
end
```

### Multi-sine Excitation

```julia
# Multiple frequencies simultaneously
load_vector = (t, dof) -> begin
    f = zeros(length(dof))
    frequencies = [0.5, 1.0, 2.0, 5.0]
    amplitudes = [1000, 800, 600, 400]
    signal = sum(A * sin(2π * freq * t) for (A, freq) in zip(amplitudes, frequencies))
    f[target_dof] = signal
    return f
end
```

### Moving Load

```julia
# Load progressing along bridge
load_vector = (t, dof) -> begin
    f = zeros(length(dof))
    bridge_progress = mod(t / crossing_time, 1.0)
    load_position = bridge_progress * bridge_length
    node = round(Int, load_position / bridge_length * n_elements) + 1
    node = clamp(node, 1, n_nodes)
    y_dof = 3 * (node - 1) + 2
    f[y_dof] = load_magnitude
    return f
end
```

---

## Constants and Defaults

### Default Values

- **Damping ratio**: 2% for all modes
- **Temperature**: 20°C reference
- **Integration tolerances**: 1e-6 for eigenvalue problems
- **Animation**: 100 frames, 20 fps

### Material Properties

Typical values for structural materials:

```julia
# Steel
ρ_steel = 7800.0    # kg/m³
E_steel = 207e9     # Pa at 20°C

# Concrete  
ρ_concrete = 2400.0 # kg/m³
E_concrete = 30e9   # Pa at 20°C

# Aluminum
ρ_aluminium = 2700.0 # kg/m³  
E_aluminium = 70e9   # Pa at 20°C
```

---

**Previous**: [Theory](theory.html) | **Next**: [Examples](examples.html)
