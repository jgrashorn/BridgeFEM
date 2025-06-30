"""
    beam_modal_ode!(du, u, p, t)

Modal space ordinary differential equation for bridge dynamics with temperature-dependent properties.

This function implements the second-order modal equations of motion:
```
q̈ᵢ + 2ζᵢωᵢq̇ᵢ + ωᵢ²qᵢ = Φᵢᵀf(t)
```

# Arguments
- `du::Vector`: Derivative vector to be filled [q̇; q̈]
- `u::Vector`: State vector [q; q̇] where q are modal displacements
- `p::NamedTuple`: Parameters containing:
  - `T_func`: Temperature function T(t)
  - `n_modes`: Number of retained modes
  - `ω_interp`: Natural frequency interpolation function
  - `ζ`: Vector of damping ratios per mode
  - `Φ_interp`: Mode shape interpolation function
  - `load_vector`: External loading function f(t, dofs)
  - `n_dofs`: Total number of DOFs
- `t::Float64`: Current time

# Implementation Details
- Natural frequencies are interpolated based on current temperature
- Mode shapes are interpolated for accurate force projection
- Damping is assumed proportional (modal damping ratios)
- Compatible with DifferentialEquations.jl solvers

# See Also
- [`solve_dynamics`](@ref): High-level dynamic simulation interface
- [`decompose_matrices`](@ref): Modal decomposition preparation
"""
function beam_modal_ode!(du, u, p, t)
    T = p.T_func(t)
    q     = u[1:p.n_modes]        # modal displacements
    qdot  = u[p.n_modes+1:end]    # modal velocities

    # Interpolate natural frequencies at current T (convert Hz to rad/s)
    ω = 2π .* p.ω_interp(T)

    # Damping ratios (could be constant, interpolated, or Rayleigh-like)
    ζ = p.ζ                      # vector of damping ratios per mode

    # Interpolate mode shapes at T
    Φ = p.Φ_interp(T)

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

function beam_physical_ode!(du, u, p, t)
    # Unpack state vector
    # u_ = u[1:p.n_dofs]        # modal displacements
    # udot_ = u[p.n_dofs+1:end] # modal velocities

    # Interpolate natural frequencies and mode shapes
    M = p.M_interp(p.T_func(t))
    ζ = p.ζ
    K = p.K_interp(p.T_func(t))

    # Assemble global load vector
    @show f = p.load_vector(t, 1:p.n_dofs)

    # Construct system matrix blocks
    Z = zeros(p.n_dofs, p.n_dofs)
    I = Matrix{Float64}(LinearAlgebra.I, p.n_dofs, p.n_dofs)
    Minv = pinv(M)
    D = 2 .* ζ .* sqrt.(diag(K) ./ diag(M)) # Rayleigh-like modal damping (if needed, adjust as appropriate)
    Dmat = Diagonal(D)
    @show A = [Z I; -Minv*K -Minv*Dmat]

    b = [zeros(p.n_dofs); Minv*f]
    du = A * u + b

end