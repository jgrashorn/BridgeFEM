using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using SparseArrays
using JSON
using Arpack
using Dates

"""
Dictionary mapping boundary condition type strings to DOF indices.

# Available constraint types:
- `"all"`: [1,2,3] - All DOFs (x, y, rotation)
- `"trans"`: [1,2] - Both translational DOFs
- `"x"`: 1 - X-direction translational DOF
- `"y"`: 2 - Y-direction translational DOF
- `"ϕ"`: 3 - Rotational DOF
"""
BCTypes = Dict(
    "all" => [1,2,3], # All DOFs (x, y, rotation)
    "trans" => [1,2], # Both translational DOFs
    "x" => 1,    # X-direction translational DOF
    "y" => 2,    # Y-direction translational DOF
    "ϕ" => 3,    # Rotational DOF
)

"""
    BridgeBC(conditions)

Boundary condition specification for bridge nodes.

# Arguments
- `conditions::Vector{Vector{Any}}`: List of boundary conditions, where each condition is 
  `[node_number, constraint_type]`

# Constraint types
Can use either string identifiers or explicit DOF vectors:
- String: `"all"`, `"trans"`, `"x"`, `"y"`, `"ϕ"` (see `BCTypes`)
- Vector: `[1, 2, 3]` for explicit DOF specification

# Examples
```julia
# Using string identifiers
bc = BridgeBC([
    [1, "all"],      # Fix all DOFs at node 1
    [10, "trans"],   # Fix translations at node 10
    [20, "y"]        # Fix y-translation at node 20
])

# Using explicit DOF vectors
bc = BridgeBC([
    [1, [1,2,3]],    # Fix all DOFs at node 1
    [5, [2]]         # Fix y-translation at node 5
])
```
"""
struct BridgeBC
    conds::Vector{Vector{Any}}

    function BridgeBC(conds::Vector{Vector{Any}})
        # If the second entry is a String, map it; otherwise, use as is
        c_ = [
            (isa(c[2], String) ? [c[1], BCTypes[c[2]]] : c)
            for c in conds
        ]
        return new(c_)
    end
end

"""
    BridgeOptions(n_elem, bc_nodes, L, ρ, A, I, E_T, cutoff_freq)

Main structure defining a bridge finite element model.

# Arguments
- `n_elem::Int`: Number of finite elements
- `bc_nodes::BridgeBC`: Boundary conditions
- `L::Float64`: Bridge length (m)
- `ρ::Float64`: Material density (kg/m³)
- `A::Float64`: Cross-sectional area (m²)
- `I::Float64`: Moment of inertia (m⁴)
- `E_T::Matrix{Float64}`: Temperature-Young's modulus data [T E; ...]
- `cutoff_freq::Float64`: Maximum frequency for modal analysis (Hz)

# Fields
- `n_nodes::Int`: Number of nodes (computed as n_elem + 1)
- `n_dofs::Int`: Number of degrees of freedom (3 per node)
- `E::Function`: Young's modulus interpolation function E(T)

# Examples
```julia
# Temperature-dependent steel bridge
E_data = [
    -10.0  250e9;   # E at -10°C
     20.0  207e9;   # E at 20°C  
     50.0  150e9    # E at 50°C
]

bc = BridgeBC([[1, "all"], [51, "y"]])
bridge = BridgeOptions(50, bc, 300.0, 7800.0, 4.0, 3.0, E_data, 50.0)
```
"""
mutable struct BridgeOptions
    n_elem::Int # Number of finite elements
    n_nodes::Int # Number of nodes
    n_dofs::Int # Number of degrees of freedom (3 per node)
    bc_nodes::BridgeBC # Boundary condition nodes (default: empty)
    L::Float64 # Beam length (m)
    ρ::Float64 # Density (kg/m^3)
    A::Float64 # Cross-section area (m^2)
    I::Float64 # Moment of inertia (m^4)
    E_T::Matrix{Float64} # Young's modulus at a specific temperature (Pa)
    E::Function # Young's modulus as a function of temperature (Pa)
    cutoff_freq::Float64 # Cutoff frequency for modes (Hz)

    function BridgeOptions(n_elem::Int, bc_nodes::BridgeBC, L::Float64, ρ::Float64, A::Float64, I::Float64, E_T::Matrix{Float64}, cutoff_freq::Float64)
        E_interp = interpolate((E_T[:,1],), E_T[:,2], Gridded(Linear()))
        E = T -> E_interp(T)
        return BridgeOptions(n_elem, bc_nodes, L, ρ, A, I, E_T, E, cutoff_freq)
    end

    function BridgeOptions(n_elem::Int, bc_nodes::BridgeBC, L::Float64, ρ::Float64, A::Float64, I::Float64, E_T::Matrix{Float64}, E::Function, cutoff_freq::Float64)
        n_nodes  = n_elem + 1
        n_dofs   = 3 * n_nodes  # 3 DOFs per node (u, v, theta)
        return new(n_elem, n_nodes, n_dofs, bc_nodes, L, ρ, A, I, E_T, E, cutoff_freq)
    end
end

function bridge_options_to_dict(bo::BridgeOptions)
    return Dict(
        "n_elem" => bo.n_elem,
        "bc_nodes" => [[c[1], c[2]] for c in bo.bc_nodes.conds],
        "L" => bo.L,
        "ρ" => bo.ρ,
        "A" => bo.A,
        "I" => bo.I,
        "E_T" => [bo.E_T[i, :] for i in 1:size(bo.E_T, 1)],
        "cutoff_freq" => bo.cutoff_freq
    )
end

function dict_to_bridge_options(dict::Dict)
    E_T_mat = Float64.(reduce(vcat, [row' for row in dict["E_T"]]))
    bconds = BridgeBC([c for c in dict["bc_nodes"]])
    return BridgeOptions(
        dict["n_elem"],
        bconds,
        dict["L"],
        dict["ρ"],
        dict["A"],
        dict["I"],
        E_T_mat,
        dict["cutoff_freq"],
    )
end

"""
    SupportElement(connection_node, connection_dofs, angle, n_elem, A, I, E_T, L, bc_bottom)

Auxiliary support structure (pier, cable, etc.) connected to the main bridge.

# Arguments
- `connection_node::Int`: Bridge node to connect to
- `connection_dofs::Vector{Int}`: DOFs to connect [1,2,3] = [x,y,θ]
- `angle::Float64`: Orientation angle in degrees (0° = horizontal right, -90° = vertical down)
- `n_elem::Int`: Number of elements in support
- `A::Float64`: Cross-sectional area (m²)
- `I::Float64`: Moment of inertia (m⁴)
- `E_T::Matrix{Float64}`: Temperature-dependent Young's modulus data
- `L::Float64`: Support length (m)
- `bc_bottom::Vector{Int}`: Constrained DOFs at support base

# Coordinate System
- **Local coordinates**: Support extends from (0,0) to (L,0) horizontally
- **Global coordinates**: Rotated by `angle` and connected to bridge
- **Connection**: First node connects to bridge, last node is constrained

# Examples
```julia
# Vertical pier with temperature dependence
pier = SupportElement(
    26,              # Connect to node 26
    [1, 2, 3],       # Connect all DOFs  
    -90.0,           # Vertical downward
    5,               # 5 elements
    0.5,             # Cross-sectional area
    0.02,            # Moment of inertia
    E_data,          # Same temperature dependence as bridge
    50.0,            # 50m height
    [1, 2, 3]        # Fix all DOFs at base
)

# Pinned support (forces only, no moments)
pin = SupportElement(5, [1, 2], -90.0, 3, 0.2, 0.01, E_data, 20.0, [1, 2])
```
"""
mutable struct SupportElement
    connection_node::Int        # Node on main bridge this support connects to
    connection_dofs::Vector{Int} # Specific DOFs to connect (e.g., [1,2] for x,y only)
    angle::Float64             # Local to global rotation angle in degrees
    n_elem::Int               # Number of elements in the support
    A::Float64                # Cross-sectional area (constant)
    I::Float64                # Moment of inertia (constant)
    E_T::Matrix{Float64}      # Young's modulus vs temperature data
    E::Function               # Young's modulus as function of temperature
    L::Float64
    bc_bottom::Vector{Int}    # Boundary conditions at bottom of support (DOF types)
end

# Constructor with temperature-dependent Young's modulus
function SupportElement(connection_node::Int, connection_dofs::Vector{Int}, angle::Float64, 
                       n_elem::Int, A::Float64, I::Float64, E_T::Matrix{Float64}, L::Float64, 
                       bc_bottom::Vector{Int})
    E_interp = interpolate((E_T[:,1],), E_T[:,2], Gridded(Linear()))
    E = T -> E_interp(T)
    return SupportElement(connection_node, connection_dofs, angle, n_elem, A, I, E_T, E, L, bc_bottom)
end

# Constructor with constant Young's modulus (creates temperature-independent support)
function SupportElement(connection_node::Int, connection_dofs::Vector{Int}, angle::Float64, 
                       n_elem::Int, E_const::Float64, A::Float64, I::Float64, L::Float64, 
                       bc_bottom::Vector{Int})
    # Create dummy temperature data for constant E
    E_T = [-100.0 E_const; 100.0 E_const]

    return SupportElement(connection_node, connection_dofs, angle, n_elem, A, I, E_T, L, bc_bottom)
end

struct SimulationOptions
    bridge::BridgeOptions
    supports::Vector{SupportElement}
    temperatures::Vector{Float64}
    damping_ratio::Float64
    total_dofs::Int
    total_elements::Int
    support_dof_mapping::Vector{Vector{Int}} # Maps support DOFs to global DOFs
    created_at::String
end

# Constructor for creating options before simulation
function SimulationOptions(bridge::BridgeOptions, supports::Vector{SupportElement}, temperatures::Vector{Float64}; damping_ratio=0.02)
    # Compute system properties
    support_dof_mapping, total_dofs = create_support_dof_mapping(bridge, supports)
    total_elements = bridge.n_elem + sum(s.n_elem for s in supports)
    
    return SimulationOptions(
        bridge, supports, temperatures, damping_ratio, total_dofs, total_elements, support_dof_mapping, string(now())
    )
end

function SimulationOptions(bridge::BridgeOptions, supports::Vector{SupportElement}, temperatures::Vector{Float64}, damping_ratio::Float64, total_dofs::Int, total_elements::Int, created_at::String)

    return SimulationOptions(bridge, supports, temperatures, damping_ratio, total_dofs, total_elements, support_dof_mapping, created_at)

end

function simulation_options_to_dict(opts::SimulationOptions)
    return Dict(
        "bridge" => bridge_options_to_dict(opts.bridge),
        "supports" => [support_element_to_dict(s) for s in opts.supports],
        "temperatures" => opts.temperatures,
        "damping_ratio" => opts.damping_ratio,
        "total_dofs" => opts.total_dofs,
        "total_elements" => opts.total_elements,
        "created_at" => opts.created_at,
    )
end

function support_element_to_dict(se::SupportElement)
    return Dict(
        "connection_node" => se.connection_node,
        "connection_dofs" => se.connection_dofs,
        "angle" => se.angle,
        "n_elem" => se.n_elem,
        "A" => se.A,
        "I" => se.I,
        "E_T" => [se.E_T[i, :] for i in 1:size(se.E_T, 1)],  # Save temperature-E data
        "L" => se.L,
        "bc_bottom" => se.bc_bottom
    )
end

function dict_to_support_element(dict::Dict)
    E_T_mat = Float64.(reduce(vcat, [row' for row in dict["E_T"]]))
    return SupportElement(
        dict["connection_node"],
        Vector{Int}(dict["connection_dofs"]),
        dict["angle"],
        dict["n_elem"],
        dict["A"],
        dict["I"],
        E_T_mat,  # Use temperature-dependent data
        dict["L"],
        Vector{Int}(dict["bc_bottom"])
    )
end

function load_simulation_options(filename::String)::SimulationOptions
    dict_data = JSON.parsefile(filename)
    
    # Reconstruct bridge options
    bridge = dict_to_bridge_options(dict_data["bridge"])
    
    # Reconstruct support elements
    supports = [dict_to_support_element(s) for s in dict_data["supports"]]
    
    # Extract other parameters
    temperatures = Vector{Float64}(dict_data["temperatures"])
    
    return SimulationOptions(
        bridge,
        supports,
        temperatures,
        dict_data["damping_ratio"],
        dict_data["total_dofs"],
        dict_data["total_elements"],
        dict_data["created_at"],
    )
end

function save_simulation_options(opts::SimulationOptions, filename::String)
    dict_data = simulation_options_to_dict(opts)
    open(filename, "w") do f
        JSON.print(f, dict_data, 2)  # Pretty print with 2-space indentation
    end
    @info "Simulation options saved to: $filename"
end

"""
    frame_elem_stiffness(EA, EI, L_e)

Compute element stiffness matrix for 2D frame element.

Returns the 6×6 local stiffness matrix for a frame element with 3 DOFs per node:
- DOFs 1,4: Axial (u₁, u₂)
- DOFs 2,5: Transverse (v₁, v₂)  
- DOFs 3,6: Rotation (θ₁, θ₂)

# Arguments
- `EA::Float64`: Axial stiffness (Young's modulus × area)
- `EI::Float64`: Flexural stiffness (Young's modulus × moment of inertia)
- `L_e::Float64`: Element length

# Returns
- `k::Matrix{Float64}`: 6×6 element stiffness matrix

# Stiffness Matrix Form
```
k = [EA/L    0      0     -EA/L    0      0   ]
    [0    12EI/L³  6EI/L²    0   -12EI/L³ 6EI/L²]
    [0     6EI/L²  4EI/L     0    -6EI/L² 2EI/L ]
    [-EA/L   0      0      EA/L    0      0   ]
    [0   -12EI/L³ -6EI/L²    0    12EI/L³ -6EI/L²]
    [0     6EI/L²  2EI/L     0    -6EI/L² 4EI/L ]
```

# Theory
Based on Euler-Bernoulli beam theory with:
- Linear axial deformation: u(x) = N₁u₁ + N₂u₂
- Cubic transverse deformation: v(x) = H₁v₁ + H₂θ₁ + H₃v₂ + H₄θ₂

# See Also
- [`frame_elem_mass`](@ref): Corresponding mass matrix
- [`assemble_stiffness!`](@ref): Global assembly process
"""
function frame_elem_stiffness(EA, EI, L_e)
    c1 = EA / L_e
    c2 = 12.0 * EI / L_e^3
    c3 = 6.0 * EI / L_e^2
    c4 = 4.0 * EI / L_e
    c5 = 2.0 * EI / L_e
    
    k = [
        c1    0.0   0.0     -c1   0.0   0.0;
        0.0   c2    c3      0.0   -c2   c3;
        0.0   c3    c4      0.0   -c3   c5;
        -c1   0.0   0.0     c1    0.0   0.0;
        0.0   -c2   -c3     0.0   c2    -c3;
        0.0   c3    c5      0.0   -c3   c4
    ]
    return k
end

"""
    frame_elem_mass(ρ, A, L)

Compute consistent mass matrix for 2D frame element.

Returns the 6×6 consistent mass matrix based on the same shape functions used 
for stiffness. Includes both translational and rotational inertia effects.

# Arguments
- `ρ::Float64`: Material density (kg/m³)
- `A::Float64`: Cross-sectional area (m²)
- `L::Float64`: Element length (m)

# Returns  
- `M::Matrix{Float64}`: 6×6 element mass matrix

# Mass Matrix Form
```
M = (ρAL/420) × [140  0   0   70   0    0  ]
                [0   156 22L  0   54  -13L]
                [0   22L 4L²  0   13L -3L²]
                [70   0   0  140   0    0  ]
                [0   54  13L  0  156  -22L]
                [0  -13L -3L² 0  -22L  4L²]
```

# Theory
- Consistent mass from Hermite interpolation functions
- Couples translational and rotational DOFs
- More accurate than lumped mass for dynamic analysis
- Preserves total mass: ∫ M dx = ρAL

# Notes
- Returns dense matrix (all entries generally non-zero)
- Rotational inertia terms scaled by L² factors
- Coupling terms (22L, 13L, etc.) ensure consistency

# See Also
- [`frame_elem_stiffness`](@ref): Corresponding stiffness matrix
- [`assemble_matrices`](@ref): Global mass assembly
"""
function frame_elem_mass(ρ::Float64, A::Float64, L::Float64)
    mass_factor = ρ * A * L / 420.0
    
    M = mass_factor * [
        140   0      0     70    0      0    ;
        0     156    22L   0     54    -13L  ;
        0     22L    4L^2  0     13L   -3L^2 ;
        70    0      0     140   0      0    ;
        0     54     13L   0     156   -22L  ;
        0    -13L   -3L^2  0    -22L    4L^2
    ]
    
    return M
end

function transformation_matrix(θ)
    c = cosd(θ)
    s = sind(θ)
    T = Matrix{Float64}(LinearAlgebra.I, 6, 6)
    T[1,1] =  c; T[1,2] =  s
    T[2,1] = -s; T[2,2] =  c
    T[4,4] =  c; T[4,5] =  s
    T[5,4] = -s; T[5,5] =  c
    return T
end

function rotated_stiffness(support::SupportElement)
    ke_local = frame_elem_stiffness(support.EA, support.EI, support.L)
    Te = transformation_matrix(support.angle)
    ke_global = Te' * ke_local * Te
    return ke_global, support.dofs
end

# Assemble global stiffness matrix
function assemble_stiffness!(K, bo::BridgeOptions, EA, EI)
    dx = bo.L / bo.n_elem
    for e = 1:bo.n_elem
        ke = frame_elem_stiffness(EA, EI, dx)
        # DOFs for element e: nodes e and e+1, each with 3 DOFs
        dofs = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*e+1, 3*e+2, 3*e+3]
        K[dofs, dofs] .+= ke
    end
    return K
end

"""
    assemble_matrices(bridge, T=20.0)

Assemble global mass and stiffness matrices for bridge at given temperature.

Constructs the finite element system matrices including temperature-dependent 
material properties and applies boundary conditions.

# Arguments
- `bridge::BridgeOptions`: Bridge model parameters
- `T::Float64=20.0`: Temperature in °C for material property evaluation

# Returns
- `M::SparseMatrixCSC{Float64}`: Global mass matrix [n_dofs × n_dofs]
- `K::SparseMatrixCSC{Float64}`: Global stiffness matrix [n_dofs × n_dofs]

# Assembly Process
1. **Initialize**: Create sparse zero matrices of size n_dofs × n_dofs
2. **Mass assembly**: Lumped mass approximation per node
   - Translational: m = ρA(dx/2) for each node
   - Rotational: I = ρA(dx³/24) for slender beam assumption
3. **Stiffness assembly**: Loop over elements and assemble element matrices
   - Material stiffness: EA(T) = E(T) × A, EI(T) = E(T) × I
   - Element connectivity: nodes i and i+1 → DOFs [3i-2:3i, 3i+1:3i+3]
4. **Boundary conditions**: Apply kinematic constraints from `bridge.bc_nodes`

# Temperature Effects
- Young's modulus E(T) evaluated via interpolation function
- Both axial (EA) and flexural (EI) stiffness scale with temperature
- Mass matrix unaffected by temperature

# Boundary Condition Application
For each constrained DOF:
- Set corresponding row/column in K to zero
- Set diagonal term to 1.0
- Zero out mass matrix diagonal term

# Example
```julia
bridge = BridgeOptions(50, bc, 300.0, 7800.0, 4.0, 3.0, E_data, 50.0)
M, K = assemble_matrices(bridge, 25.0)  # At 25°C

# Check system properties
println("System size: \$(size(K))")
println("Condition number: \$(cond(Array(K)))")
```

# See Also
- [`assemble_matrices_with_supports`](@ref): Extended version with support elements
- [`apply_bc!`](@ref): Boundary condition enforcement
- [`frame_elem_stiffness`](@ref), [`frame_elem_mass`](@ref): Element matrices
"""
function assemble_matrices(bo::BridgeOptions, T::Float64=20.0)

    # Discretization
    dx = bo.L / bo.n_elem
    n_dof = bo.n_dofs    # 3 DOFs per node (u, v, theta)

    M = spzeros(n_dof, n_dof)
    K = spzeros(n_dof, n_dof)

    # Mass matrix (lumped for simplicity: translational and rotational inertia per node)
    m_trans = bo.ρ * bo.A * dx / 2  # each node shares half element mass
    m_rot   = bo.ρ * bo.A * dx^3 / 24  # rotational inertia for slender beam

    for i in 1:bo.n_nodes
        # Translational DOFs (u and v)
        M[3*(i-1)+1, 3*(i-1)+1] = m_trans  # u direction
        M[3*(i-1)+2, 3*(i-1)+2] = m_trans  # v direction
        # Rotational DOF (theta)
        M[3*(i-1)+3, 3*(i-1)+3] = m_rot    # rotation
    end

    E = bo.E(T)  # Young's modulus at temperature T
    assemble_stiffness!(K, bo, E * bo.A, E * bo.I)

    apply_bc!(M, K, bo.bc_nodes.conds)

    return M, K
end

function apply_bc!(M, K, bc_dofs)
    n_dofs = size(K, 1)
    asdf = 1
    for bc in bc_dofs
        node = bc[1]
        dof_types = bc[2]  # DOF type(s)
        dof_indices = 3 * (node - 1) .+ dof_types  # Convert to global DOF indices
        
        for d_ in dof_indices
            # @info "Applying boundary condition at DOF $d_"
            if d_ <= n_dofs  # Check bounds
                K[:, d_] .= 0.0
                K[d_, :] .= 0.0
                K[d_, d_] = 1.0
                M[d_, d_] = 0.0
            end
        end
    end
end

function assemble_local_support(support::SupportElement, T::Float64=20.0)
    n_nodes = support.n_elem + 1
    n_dofs = 3 * n_nodes
    dx = support.L / support.n_elem
    
    K_local = zeros(n_dofs, n_dofs)
    
    # Get temperature-dependent Young's modulus
    E = support.E(T)
    EA = E * support.A
    EI = E * support.I
    
    # Assemble support elements
    for e = 1:support.n_elem
        ke = frame_elem_stiffness(EA, EI, dx)
        # DOFs for element e: nodes e and e+1
        dofs = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*e+1, 3*e+2, 3*e+3]
        K_local[dofs, dofs] .+= ke
    end
    
    return K_local
end

function create_support_dof_mapping(bo::BridgeOptions, supports::Vector{SupportElement})
    global_dof_offset = bo.n_dofs
    support_dof_maps = Vector{Vector{Int}}()
    
    for (i, support) in enumerate(supports)
        n_support_nodes = support.n_elem + 1
        n_support_dofs = 3 * n_support_nodes
        
        # Initialize mapping vector
        local_to_global = zeros(Int, n_support_dofs)
        
        # SIMPLIFIED: First node (index 1) connects to bridge
        bridge_connection_node = support.connection_node
        bridge_connection_dofs = 3 * (bridge_connection_node - 1) .+ support.connection_dofs
        
        # First node of support connects to bridge
        connection_node_dofs = [1, 2, 3]  # First node DOFs
        
        # Map connected DOFs from first node to bridge
        for (j, connection_dof) in enumerate(support.connection_dofs)
            local_support_dof = connection_node_dofs[connection_dof]
            local_to_global[local_support_dof] = bridge_connection_dofs[j]
        end
        
        # Assign new DOFs to all non-connected DOFs
        for dof_idx in 1:n_support_dofs
            if local_to_global[dof_idx] == 0  # Not yet assigned
                global_dof_offset += 1
                local_to_global[dof_idx] = global_dof_offset
            end
        end
        
        push!(support_dof_maps, local_to_global)
    end
    
    total_dofs = global_dof_offset
    return support_dof_maps, total_dofs
end

function get_node_from_dof(so::SimulationOptions,dof::Int)

end

"""
    assemble_matrices_with_supports(bridge, supports, T=20.0)

Assemble expanded system matrices including bridge and support structures.

# Arguments
- `bridge::BridgeOptions`: Main bridge configuration
- `supports::Vector{SupportElement}`: Support structures
- `T::Float64`: Temperature (°C) for material properties

# Returns
- `M::SparseMatrixCSC`: Expanded global mass matrix
- `K::SparseMatrixCSC`: Expanded global stiffness matrix

# Method
1. Creates DOF mapping between bridge and support systems
2. Assembles bridge matrices in upper-left block
3. For each support:
   - Assembles local support matrices
   - Rotates to global coordinates using transformation matrix
   - Maps to expanded global DOF numbering
   - Applies boundary conditions at support base
"""
function assemble_matrices_with_supports(bo::BridgeOptions, supports::Vector{SupportElement}, T::Float64=20.0)
    # Get DOF mappings
    support_dof_maps, total_dofs = create_support_dof_mapping(bo, supports)
    
    # Initialize expanded matrices
    M = spzeros(total_dofs, total_dofs)
    K = spzeros(total_dofs, total_dofs)
    
    # Assemble main bridge (already temperature-dependent)
    M_bridge, K_bridge = assemble_matrices(bo, T)
    M[1:bo.n_dofs, 1:bo.n_dofs] .= M_bridge
    K[1:bo.n_dofs, 1:bo.n_dofs] .= K_bridge
    
    # Assemble each support (now temperature-dependent)
    for (i, support) in enumerate(supports)
        # Get local support matrices at temperature T
        K_local = assemble_local_support(support, T)  # FIXED: Pass temperature
        M_local = create_support_mass_matrix(support, bo.ρ)  # Mass not temperature dependent
        
        # Rotate BOTH matrices to global coordinates
        n_support_nodes = support.n_elem + 1
        T_expanded = create_expanded_transformation(support.angle, n_support_nodes)
        
        K_rotated = T_expanded' * K_local * T_expanded
        M_rotated = T_expanded' * M_local * T_expanded
        
        # Map to global DOFs
        dof_map = support_dof_maps[i]
        K[dof_map, dof_map] .+= K_rotated
        M[dof_map, dof_map] .+= M_rotated
        
        # Apply boundary conditions to LAST node (fixed base)
        n_support_nodes = support.n_elem + 1
        fixed_node_local_dofs = 3*(n_support_nodes-1) .+ [1, 2, 3]  # Last node
        fixed_node_global_dofs = dof_map[fixed_node_local_dofs]
        
        # Apply boundary conditions
        for dof_type in support.bc_bottom
            global_dof = fixed_node_global_dofs[dof_type]
            K[:, global_dof] .= 0.0
            K[global_dof, :] .= 0.0
            K[global_dof, global_dof] = 1.0
            M[:, global_dof] .= 0.0
            M[global_dof, :] .= 0.0
            M[global_dof, global_dof] = 0.0
        end
    end
    
    return M, K
end

function create_expanded_transformation(angle::Float64, n_nodes::Int)
    n_dofs = 3 * n_nodes
    T_expanded = Matrix{Float64}(LinearAlgebra.I, n_dofs, n_dofs)
    
    T_single = transformation_matrix(angle)
    
    for node = 1:n_nodes
        dof_start = 3 * (node - 1) + 1
        dof_end = dof_start + 2
        
        # Apply 2D rotation to x,y DOFs (rotation DOF unchanged)
        T_expanded[dof_start:dof_start+1, dof_start:dof_start+1] = T_single[1:2, 1:2]
    end
    
    return T_expanded
end

function create_support_mass_matrix(support::SupportElement, ρ::Float64)
    n_nodes = support.n_elem + 1
    n_dofs = 3 * n_nodes
    dx = support.L / support.n_elem
    
    M_local = zeros(n_dofs, n_dofs)
    
    # Mass matrix is not temperature dependent (density and geometry constant)
    # Use reference area for mass calculation
    A_ref = support.A
    
    # Assemble mass matrix for each element
    for e = 1:support.n_elem
        me = frame_elem_mass(ρ, A_ref, dx)
        dofs = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*e+1, 3*e+2, 3*e+3]
        M_local[dofs, dofs] .+= me
    end
    
    return M_local
end