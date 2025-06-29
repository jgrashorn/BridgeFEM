using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using SparseArrays
using JSON
using Arpack
using Dates

BCTypes = Dict(
    "all" => [1,2,3], # All DOFs (x, y, rotation)
    "trans" => [1,2], # Both translational DOFs
    "x" => 1,    # X-direction translational DOF
    "y" => 2,    # Y-direction translational DOF
    "ϕ" => 3,    # Rotational DOF
)

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

    function BridgeOptions(n_elem::Int, bc_nodes::BridgeBC, L::Float64, ρ::Float64, A::Float64, I::Float64, E::Function, T::Vector{Float64}, cutoff_freq::Float64)
        E_T = hcat(T, E.(T))  # Create a matrix of temperatures and corresponding Young's moduli
        return BridgeOptions(n_elem, bc_nodes, L, ρ, A, I, E_T, E, cutoff_freq)
    end

    function BridgeOptions(n_elem::Int, bc_nodes::BridgeBC, L::Float64, ρ::Float64, A::Float64, I::Float64, E::Function, cutoff_freq::Float64)
        @warn "No temperature provided. When saving options, please provide temperatures!"
        return BridgeOptions(n_elem, bc_nodes, L, ρ, A, I, Matrix{Float64}(undef, 0,2), E, cutoff_freq)
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

function save_options_to_json(bo::BridgeOptions; fname::String="data/bridge_options.json", T::Union{Vector{Float64},Nothing}=nothing)

    if bo.E_T == Matrix{Float64}(undef, 0,2) && !isnothing(T)
        bo.E_T = hcat(T, bo.E.(T))  # Create a matrix of temperatures and corresponding Young's moduli
    elseif bo.E_T == Matrix{Float64}(undef, 0,2) && isnothing(T)
        @error "No temperature provided. Cannot save to JSON!\n Please provide temperatures as a Vector{Float64} to the `T` keyword argument."
        return nothing
    end
    options_dict = Dict(
        "n_elem" => bo.n_elem,
        "bcconds" => [[c[1], c[2]] for c in bo.bc_nodes.conds],
        "L" => bo.L,
        "ρ" => bo.ρ,
        "A" => bo.A,
        "I" => bo.I,
        "E_T" => [bo.E_T[i, :] for i in 1:size(bo.E_T, 1)],
        "cutoff_freq" => bo.cutoff_freq
    )

    # Ensure the directory exists before saving
    dir = dirname(fname)
    if !isdir(dir)
        mkpath(dir)
    end
    open(fname, "w") do file
        JSON.print(file, options_dict, 4)
        flush(file)
    end
    @info "Bridge options saved to $fname"

end

function load_options_from_json(fname::String="data/bridge_options.json")
    options_dict = JSON.parsefile(fname)
    E_T_mat = Float64.(reduce(vcat, [row' for row in options_dict["E_T"]]))
    bconds = BridgeBC([c for c in options_dict["bcconds"]])
    return BridgeOptions(
        options_dict["n_elem"],
        bconds,
        options_dict["L"],
        options_dict["ρ"],
        options_dict["A"],
        options_dict["I"],
        E_T_mat,
        options_dict["cutoff_freq"],
    )
end

struct SupportElement
    connection_node::Int        # Node on main bridge this support connects to
    connection_dofs::Vector{Int} # Specific DOFs to connect (e.g., [1,2] for x,y only)
    angle::Float64             # Local to global rotation angle in degrees
    n_elem::Int               # Number of elements in the support
    EA::Float64
    EI::Float64
    L::Float64
    bc_bottom::Vector{Int}    # Boundary conditions at bottom of support (DOF types)
end

struct SimulationOptions
    # Bridge parameters
    bridge::BridgeOptions
    
    # Support elements
    supports::Vector{SupportElement}
    
    # Temperature range
    temperatures::Vector{Float64}
    
    # Dynamic simulation parameters
    damping_ratio::Float64
    
    # Computed properties (for reconstruction)
    total_dofs::Int
    total_elements::Int
    
    # Metadata
    created_at::String
end

# Constructor for creating options before simulation
function SimulationOptions(bridge::BridgeOptions, supports::Vector{SupportElement}, temperatures::Vector{Float64}; damping_ratio=0.02)
    # Compute system properties
    _, total_dofs = create_support_dof_mapping(bridge, supports)
    total_elements = bridge.n_elem + sum(s.n_elem for s in supports)
    
    return SimulationOptions(
        bridge, supports, temperatures, damping_ratio, total_dofs, total_elements, string(now())
    )
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
        "EA" => se.EA,
        "EI" => se.EI,
        "L" => se.L,
        "bc_bottom" => se.bc_bottom
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

function dict_to_support_element(dict::Dict)
    return SupportElement(
        dict["connection_node"],
        Vector{Int}(dict["connection_dofs"]),
        dict["angle"],
        dict["n_elem"],
        dict["EA"],
        dict["EI"],
        dict["L"],
        Vector{Int}(dict["bc_bottom"])
    )
end

function save_simulation_options(opts::SimulationOptions, filename::String)
    dict_data = simulation_options_to_dict(opts)
    open(filename, "w") do f
        JSON.print(f, dict_data, 2)  # Pretty print with 2-space indentation
    end
    @info "Simulation options saved to: $filename"
end

function rotated_stiffness(support::SupportElement)
    ke_local = frame_elem_stiffness(support.EA, support.EI, support.L)
    Te = transformation_matrix(support.angle)
    ke_global = Te' * ke_local * Te
    return ke_global, support.dofs
end

# Element stiffness matrix for 2D frame element (3 DOFs per node: u, v, theta)
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

function decompose_matrices(M, K)
    # Collect mode shapes at different temperatures
    n_dof = size(M, 1)
    n_temps = size(M, 3)
    n_mode = n_dof # Number of modes to keep

    Φ_tensor = zeros(n_dof, n_mode, n_temps)
    Φ_tensor_unnormalized = zeros(n_dof, n_mode, n_temps)
    ω_matrix = zeros(n_mode, n_temps)

    decomp = [eigen(K[:,:,i], M[:,:,i]) for i in axes(M, 3)]

    λs_all = [d.values for d in decomp]
    vecs_all = [d.vectors for d in decomp]

    # decomp = [eigs(K[:,:,i], M[:,:,i], nev=n_mode, which=:SM) for i in axes(M, 3)]
    # λs_all = [d[1] for d in decomp]
    # vecs_all = [d[2] for d in decomp]

    for i in axes(M, 3)
        λs = λs_all[i]
        vecs = vecs_all[i]
        # Normalize mode shapes (optional)
        for i in 1:n_mode
            vecs[:,i] ./= norm(vecs[:,i])
        end

        Φ_tensor[:,:,i] .= vecs

        if i > 1
            for dof in 1:n_mode
                dot(Φ_tensor[:,dof,i], Φ_tensor[:,dof,i-1])
                if dot(Φ_tensor[:,dof,i], Φ_tensor[:,dof,i-1]) < 0
                    # @info "Flipping mode shape $i"
                    Φ_tensor[:,dof,i] .*= -1
                end
            end
        end

        ωs = sqrt.(real(λs))./ (2π)  # Convert eigenvalues to natural frequencies (Hz)
        Φ_tensor_unnormalized[:,:,i] .= Φ_tensor[:,:,i]

        # Mass-normalize the mode shapes
        for j in axes(Φ_tensor, 2)
            mi = Φ_tensor[:, j, i]' * M[:,:,i] * Φ_tensor[:, j, i]
            Φ_tensor[:, j, i] ./= sqrt(mi)
        end

        ω_matrix[:,i] .= ωs
    end

    return ω_matrix, Φ_tensor, Φ_tensor_unnormalized  # Return both normalized and unnormalized mode shapes
end

function assemble_and_decompose(bo::BridgeOptions, Ts::Vector{Float64}; supports::Vector{SupportElement} = SupportElement[])
    nTs = length(Ts)
    @info "Assembling matrices with supports for $nTs temperatures"
    mats = [assemble_matrices_with_supports(bo, supports, T) for T in Ts]
    M = cat([mats[i][1] for i in 1:nTs]..., dims=3)
    K = cat([mats[i][2] for i in 1:nTs]..., dims=3)

    @info "Decomposing matrices"
    λs, vectors, vectors_unnormalized = decompose_matrices(M, K)

    keep_modes = λs[:,1] .< bo.cutoff_freq
    λs = λs[keep_modes, :]
    vectors = vectors[:, keep_modes, :]
    vectors_unnormalized = vectors_unnormalized[:, keep_modes, :]

    return M, K, λs, vectors, vectors_unnormalized
end

# function assemble_and_decompose(bo::BridgeOptions, Ts::Vector{Float64})
#     nTs = length(Ts)
#     @info "Assembling matrices for $nTs temperatures"
#     mats = [assemble_matrices(bo, T) for T in Ts]
#     M = cat([mats[i][1] for i in 1:nTs]..., dims=3)
#     K = cat([mats[i][2] for i in 1:nTs]..., dims=3)

#     @info "Decomposing matrices"
#     λs, vectors, vectors_unnormalized = decompose_matrices(M, K)

#     keep_modes = λs[:,1] .< bo.cutoff_freq
#     λs = λs[keep_modes, :]
#     vectors = vectors[:, keep_modes, :]
#     vectors_unnormalized = vectors_unnormalized[:, keep_modes, :]

#     return M, K, λs, vectors, vectors_unnormalized
# end

function beam_modal_ode!(du, u, p, t)
    T = p.T_func(t)
    q     = u[1:p.n_modes]        # modal displacements
    qdot  = u[p.n_modes+1:end]    # modal velocities

    # Interpolate natural frequencies at current T (convert Hz to rad/s)
    ω = 2π .* p.ω_interp(1:p.n_modes, T)

    # Damping ratios (could be constant, interpolated, or Rayleigh-like)
    ζ = p.ζ                      # vector of damping ratios per mode

    # Interpolate mode shapes at T
    Φ = p.Φ_interp(1:p.n_dofs, 1:p.n_modes, T)

    # Assemble global load vector
    f = p.load_vector(t,1:p.n_dofs)

    # Project force onto each mode
    fhat = Φ' * f

    # Modal accelerations
    qddot = [-2ζ[i]*ω[i]*qdot[i] - (ω[i]^2)*q[i] + fhat[i] for i in 1:p.n_modes]

    # Fill derivative vector
    du[1:p.n_modes] .= qdot
    du[p.n_modes+1:end] .= qddot
end

function reconstruct_physical(bo::BridgeOptions, q_full, Φ_interp, T_func, time)
    n_dofs  = bo.n_dofs
    n_modes_total = size(q_full, 1)
    n_modes = n_modes_total ÷ 2
    n_times = length(time)

    u_full  = zeros(n_dofs, n_times)
    du_full = zeros(n_dofs, n_times)

    for (i, t) in enumerate(time)
        T_now = T_func(t)
        Φ = Φ_interp(1:n_dofs, 1:n_modes, T_now)

        q_disp = q_full[1:n_modes, i]
        q_vel  = q_full[n_modes+1:end, i]

        u_full[:, i]  .= Φ * q_disp
        du_full[:, i] .= Φ * q_vel
    end

    return u_full, du_full
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

function assemble_local_support(support::SupportElement)
    n_nodes = support.n_elem + 1
    n_dofs = 3 * n_nodes
    dx = support.L / support.n_elem
    
    K_local = zeros(n_dofs, n_dofs)
    
    # Assemble support elements
    for e = 1:support.n_elem
        ke = frame_elem_stiffness(support.EA, support.EI, dx)
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

function assemble_matrices_with_supports(bo::BridgeOptions, supports::Vector{SupportElement}, T::Float64=20.0)
    # Get DOF mappings
    support_dof_maps, total_dofs = create_support_dof_mapping(bo, supports)
    
    # Initialize expanded matrices
    M = spzeros(total_dofs, total_dofs)
    K = spzeros(total_dofs, total_dofs)
    
    # Assemble main bridge
    M_bridge, K_bridge = assemble_matrices(bo, T)
    M[1:bo.n_dofs, 1:bo.n_dofs] .= M_bridge
    K[1:bo.n_dofs, 1:bo.n_dofs] .= K_bridge
    
    # Assemble each support
    for (i, support) in enumerate(supports)
        # Get local support matrices
        K_local = assemble_local_support(support)
        M_local = create_support_mass_matrix(support, bo.ρ)  # Add this function if missing
        
        # Rotate BOTH matrices to global coordinates
        n_support_nodes = support.n_elem + 1
        T_expanded = create_expanded_transformation(support.angle, n_support_nodes)
        
        K_rotated = T_expanded' * K_local * T_expanded
        M_rotated = T_expanded' * M_local * T_expanded  # MISSING: Rotate mass matrix too!
        
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
            # FIXED: Zero out entire row and column (not just diagonal)
            K[:, global_dof] .= 0.0    # Zero column
            K[global_dof, :] .= 0.0    # Zero row
            K[global_dof, global_dof] = 1.0  # Set diagonal
            
            M[:, global_dof] .= 0.0    # Zero column  
            M[global_dof, :] .= 0.0    # Zero row
            M[global_dof, global_dof] = 0.0  # Zero diagonal (no mass for constrained DOF)
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
    
    # Assemble mass matrix for each element
    for e = 1:support.n_elem
        # Consistent mass matrix for frame element
        me = frame_elem_mass(ρ, support.EA/support.EI, dx)  # You'll need this function
        
        # Or simpler lumped mass:
        # me = lumped_frame_mass(ρ, A_support, dx)
        
        dofs = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*e+1, 3*e+2, 3*e+3]
        M_local[dofs, dofs] .+= me
    end
    
    return M_local
end

function frame_elem_mass(ρ::Float64, A::Float64, L::Float64)
    # Consistent mass matrix for 2D frame element
    # This is a standard result from FEM textbooks
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

function plot_bridge_with_supports(bo::BridgeOptions, supports::Vector{SupportElement}; 
                                 mode_shape=nothing, scale_factor=1.0, T=20.0,
                                 title="Bridge Structure", show_nodes=true)
    
    # Bridge nodes (assuming horizontal bridge)
    bridge_x = collect(range(0, bo.L, length=bo.n_nodes))
    bridge_y = zeros(bo.n_nodes)
    
    # If mode shape is provided, add deformation DIRECTLY
    if mode_shape !== nothing
        for i in 1:bo.n_nodes
            if 3*(i-1)+1 <= length(mode_shape) && 3*(i-1)+2 <= length(mode_shape)
                bridge_x[i] += scale_factor * mode_shape[3*(i-1)+1]  # x displacement
                bridge_y[i] += scale_factor * mode_shape[3*(i-1)+2]  # y displacement
            end
        end
    end
    
    # Plot bridge
    p = plot(bridge_x, bridge_y, linewidth=3, color=:blue, label="Bridge", 
             aspect_ratio=:equal, grid=true)
    
    if show_nodes
        scatter!(bridge_x, bridge_y, color=:blue, markersize=4, label="Bridge Nodes")
    end
    
    # Plot each support
    support_dof_maps, total_dofs = create_support_dof_mapping(bo, supports)
    
    for (i, support) in enumerate(supports)
        # Connection point on bridge
        conn_node = support.connection_node
        conn_x_undeformed = bo.L * (conn_node - 1) / (bo.n_nodes - 1)  # Undeformed bridge x-position
        conn_y_undeformed = 0.0  # Undeformed bridge y-position
        
        # Support geometry
        n_support_nodes = support.n_elem + 1
        angle_rad = deg2rad(support.angle)
        
        # UNDEFORMED support positions
        support_x_undeformed = zeros(n_support_nodes)
        support_y_undeformed = zeros(n_support_nodes)
        
        for j in 1:n_support_nodes
            # Distance along support from connection point
            distance = (j-1) * support.L / (n_support_nodes - 1)
            support_x_undeformed[j] = conn_x_undeformed + distance * cos(angle_rad)
            support_y_undeformed[j] = conn_y_undeformed + distance * sin(angle_rad)
        end
        
        # DEFORMED support positions (if mode shape provided)
        support_x_deformed = copy(support_x_undeformed)
        support_y_deformed = copy(support_y_undeformed)
        
        if mode_shape !== nothing
            dof_map = support_dof_maps[i]
            
            for j in 1:n_support_nodes
                local_x_idx = 3*(j-1) + 1
                local_y_idx = 3*(j-1) + 2
                
                if local_x_idx <= length(dof_map) && local_y_idx <= length(dof_map)
                    global_x_dof = dof_map[local_x_idx]
                    global_y_dof = dof_map[local_y_idx]
                    
                    # Add deformation DIRECTLY from global mode shape
                    if global_x_dof <= length(mode_shape)
                        support_x_deformed[j] += scale_factor * mode_shape[global_x_dof]
                    end
                    if global_y_dof <= length(mode_shape)
                        support_y_deformed[j] += scale_factor * mode_shape[global_y_dof]
                    end
                end
            end
        end
        
        # Plot support
        plot!(support_x_deformed, support_y_deformed, linewidth=2, color=:red, 
              label=i==1 ? "Supports" : "")
        
        if show_nodes
            scatter!(support_x_deformed, support_y_deformed, color=:red, markersize=3, 
                    label="")
        end
        
        # Fixed base at UNDEFORMED position (last node)
        scatter!([support_x_undeformed[end]], [support_y_undeformed[end]], 
                color=:black, marker=:square, markersize=8, 
                label=i==1 ? "Fixed Base" : "")
    end
    
    xlabel!("X Position (m)")
    ylabel!("Y Position (m)")
    title!(title)
    
    return p
end

function plot_mode_shape(bo::BridgeOptions, supports::Vector{SupportElement}, 
                        mode_shape::Vector{Float64}, mode_num::Int;
                        scale_factor=50.0, T=20.0)
    
    # Plot undeformed structure
    p1 = plot_bridge_with_supports(bo, supports, T=T, 
                                  title="Undeformed Structure", show_nodes=false)
    plot!(p1, linewidth=2, alpha=0.3, color=:gray)
    
    # Plot deformed structure (mode shape)
    p2 = plot_bridge_with_supports(bo, supports, mode_shape=mode_shape, 
                                  scale_factor=scale_factor, T=T,
                                  title="Mode $mode_num (Scaled by $scale_factor)", 
                                  show_nodes=true)
    
    # Combine plots
    plot(p1, p2, layout=(1,2), size=(1200, 400))
end

function animate_mode(bo::BridgeOptions, supports::Vector{SupportElement}, 
                     mode_shape::Vector{Float64}, mode_num::Int;
                     scale_factor=50.0, T=20.0, n_frames=30, fsize = (800,800))
    
    # Calculate sensible fixed limits
    bridge_extent = bo.L
    support_extent = maximum([s.L for s in supports], init=0.0)
    max_deformation = scale_factor * maximum(abs.(mode_shape), init=0.0)
    
    padding = 0.0
    x_limits = (-padding - support_extent, bridge_extent + padding + support_extent)
    y_limits = (-padding - support_extent - max_deformation, 
                padding + max_deformation)
    
    # Create animation with fixed layout
    anim = Animation()
    
    for i in 1:n_frames
        phase = 2π * (i-1) / n_frames
        current_scale = scale_factor * sin(phase)
        
        p = plot_bridge_with_supports(bo, supports, 
                                    mode_shape=mode_shape, 
                                    scale_factor=current_scale, 
                                    T=T,
                                    title="Mode $mode_num Animation",
                                    show_nodes=true)
        
        # Apply consistent layout
        plot!(p, 
            xlims=x_limits, 
            ylims=y_limits,
            aspect_ratio=1,        # Reduce this value to make plot wider
            size=fsize,        # Increase figure width
            legend=:topright,        
            legend_background_color=:white)
        
        frame(anim)
    end
    
    return anim
end

