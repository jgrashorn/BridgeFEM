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