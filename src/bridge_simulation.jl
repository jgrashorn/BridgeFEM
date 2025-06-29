using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using SparseArrays
using JSON
using Arpack

BCTypes = Dict(
    "all" => [1,2,3], # All DOFs (x, y, rotation)
    "trans" => [1,2], # Both translational DOFs
    "x" => 1,    # X-direction translational DOF
    "y" => 2,    # Y-direction translational DOF
    "ϕ" => 3,    # Rotational DOF
    "both" => [1,3], # Legacy: x and rotation (for backward compatibility)
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
    for bc in bc_dofs
        node = bc[1]
        dof_types = bc[2]  # DOF type(s)
        dof_indices = 3 * (node - 1) .+ dof_types  # Convert to global DOF indices
        
        for d_ in dof_indices
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

        # Mass-normalize the mode shapes
        for j in axes(Φ_tensor, 2)
            mi = Φ_tensor[:, j, i]' * M[:,:,i] * Φ_tensor[:, j, i]
            Φ_tensor[:, j, i] ./= sqrt(mi)
        end

        ω_matrix[:,i] .= ωs
    end

    return ω_matrix, Φ_tensor
end

function assemble_and_decompose(bo::BridgeOptions, Ts::Vector{Float64})
    nTs = length(Ts)
    @info "Assembling matrices for $nTs temperatures"
    mats = [assemble_matrices(bo, T) for T in Ts]
    M = cat([mats[i][1] for i in 1:nTs]..., dims=3)
    K = cat([mats[i][2] for i in 1:nTs]..., dims=3)

    @info "Decomposing matrices"
    λs, vectors = decompose_matrices(M, K)

    keep_modes = λs[:,1] .< bo.cutoff_freq
    λs = λs[keep_modes, :]
    vectors = vectors[:, keep_modes, :]

    return M, K, λs, vectors
end

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