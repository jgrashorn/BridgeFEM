---
layout: default
title: Theory - BridgeFEM.jl
---

# Theoretical Background

This page describes the finite element formulation and model order reduction techniques used in BridgeFEM.jl.

## Table of Contents
- [Finite Element Method](#finite-element-method)
- [2D Frame Elements](#2d-frame-elements)
- [Temperature Dependence](#temperature-dependence)
- [Support Modeling](#support-modeling)
- [Modal Analysis](#modal-analysis)
- [Model Order Reduction](#model-order-reduction)
- [Dynamic Simulation](#dynamic-simulation)

---

## Finite Element Method

### Governing Equations

For a structural system, the equations of motion in physical coordinates are:

$$\mathbf{M}\ddot{\mathbf{u}} + \mathbf{C}\dot{\mathbf{u}} + \mathbf{K}(\mathbf{T})\mathbf{u} = \mathbf{f}(t)$$

where:
- $\mathbf{M}$ = mass matrix
- $\mathbf{C}$ = damping matrix  
- $\mathbf{K}(\mathbf{T})$ = temperature-dependent stiffness matrix
- $\mathbf{u}$ = displacement vector
- $\mathbf{f}(t)$ = external force vector

### Discretization

The bridge is discretized using 2D frame elements with 3 degrees of freedom per node:
- $u$ = horizontal displacement
- $v$ = vertical displacement  
- $\theta$ = rotation about z-axis

For $n$ elements with $n+1$ nodes, the total system size is $3(n+1)$ DOFs.

---

## 2D Frame Elements

### Element Stiffness Matrix

Each frame element has 6 DOFs (3 per node). The local element stiffness matrix is:

$$\mathbf{k}^e = \begin{bmatrix}
\frac{EA}{L} & 0 & 0 & -\frac{EA}{L} & 0 & 0 \\
0 & \frac{12EI}{L^3} & \frac{6EI}{L^2} & 0 & -\frac{12EI}{L^3} & \frac{6EI}{L^2} \\
0 & \frac{6EI}{L^2} & \frac{4EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{2EI}{L} \\
-\frac{EA}{L} & 0 & 0 & \frac{EA}{L} & 0 & 0 \\
0 & -\frac{12EI}{L^3} & -\frac{6EI}{L^2} & 0 & \frac{12EI}{L^3} & -\frac{6EI}{L^2} \\
0 & \frac{6EI}{L^2} & \frac{2EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{4EI}{L}
\end{bmatrix}$$

where:
- $E$ = Young's modulus
- $A$ = cross-sectional area
- $I$ = moment of inertia
- $L$ = element length

### Element Mass Matrix

The consistent mass matrix for a frame element is:

$$\mathbf{m}^e = \frac{\rho A L}{420} \begin{bmatrix}
140 & 0 & 0 & 70 & 0 & 0 \\
0 & 156 & 22L & 0 & 54 & -13L \\
0 & 22L & 4L^2 & 0 & 13L & -3L^2 \\
70 & 0 & 0 & 140 & 0 & 0 \\
0 & 54 & 13L & 0 & 156 & -22L \\
0 & -13L & -3L^2 & 0 & -22L & 4L^2
\end{bmatrix}$$

### Assembly Process

Global matrices are assembled using the standard finite element procedure:

```julia
function assemble_matrices(bridge, T)
    M = spzeros(n_dofs, n_dofs)
    K = spzeros(n_dofs, n_dofs)
    
    for e = 1:n_elements
        # Element matrices
        ke = frame_elem_stiffness(E(T)*A, E(T)*I, L_element)
        me = frame_elem_mass(ρ, A, L_element)
        
        # DOF connectivity
        dofs = element_dofs(e)
        
        # Assemble
        K[dofs, dofs] += ke
        M[dofs, dofs] += me
    end
    
    return M, K
end
```

---

## Temperature Dependence

### Material Property Interpolation

Young's modulus varies with temperature according to user-defined data points:

$$E(T) = \text{interpolate}(T; \{T_i, E_i\}_{i=1}^n)$$

Linear interpolation is used between data points:

```julia
E_interp = interpolate((T_data[:,1],), T_data[:,2], Gridded(Linear()))
E = T -> E_interp(T)
```

### Stiffness Matrix Update

At each temperature $T$, the stiffness matrix becomes:

$$\mathbf{K}(T) = \sum_{e=1}^{n_{elem}} \mathbf{L}_e^T \mathbf{k}^e(E(T)) \mathbf{L}_e$$

where $\mathbf{L}_e$ is the Boolean connectivity matrix for element $e$.

---

## Support Modeling

### Coordinate Transformation

Support elements can be oriented at arbitrary angles. The transformation from local to global coordinates is:

$$\mathbf{k}_{global}^{support} = \mathbf{T}^T \mathbf{k}_{local}^{support} \mathbf{T}$$

where the transformation matrix for a 2D rotation by angle $\theta$ is:

$$\mathbf{T} = \begin{bmatrix}
\cos\theta & \sin\theta & 0 & 0 & 0 & 0 \\
-\sin\theta & \cos\theta & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & \cos\theta & \sin\theta & 0 \\
0 & 0 & 0 & -\sin\theta & \cos\theta & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{bmatrix}$$

### DOF Connectivity

Support elements connect to the main bridge through shared DOFs:

1. **Bridge DOFs**: $\{u_i, v_i, \theta_i\}$ for connection node $i$
2. **Support DOFs**: Local support numbering mapped to global system
3. **Connectivity**: User-specified which DOFs are shared (e.g., $[1,2]$ for pinned, $[1,2,3]$ for rigid)

### Expanded System Assembly

The final system includes both bridge and support DOFs:

$$\begin{bmatrix}
\mathbf{M}_{bridge} & \mathbf{0} \\
\mathbf{0} & \mathbf{M}_{support}
\end{bmatrix} \ddot{\mathbf{u}} + 
\begin{bmatrix}
\mathbf{K}_{bridge} & \mathbf{K}_{coupling} \\
\mathbf{K}_{coupling}^T & \mathbf{K}_{support}
\end{bmatrix} \mathbf{u} = \mathbf{f}$$

---

## Modal Analysis

### Generalized Eigenvalue Problem

Natural frequencies and mode shapes are found by solving:

$$\mathbf{K}\boldsymbol{\phi}_i = \lambda_i \mathbf{M}\boldsymbol{\phi}_i$$

where:
- $\lambda_i = \omega_i^2$ = squared natural frequency
- $\boldsymbol{\phi}_i$ = mode shape vector
- $\omega_i = 2\pi f_i$ = angular frequency

### Solution Method

For large sparse systems, ARPACK is used for efficient eigenvalue computation:

```julia
λ, φ = eigs(K, M, nev=n_modes, which=:SM, sigma=1e-6)
frequencies = sqrt.(λ) / (2π)  # Convert to Hz
```

### Mode Normalization

Mode shapes are typically normalized such that:

$$\boldsymbol{\phi}_i^T \mathbf{M} \boldsymbol{\phi}_i = 1$$

This ensures unit modal mass for each mode.

---

## Model Order Reduction

### Modal Truncation

The physical displacement vector is approximated using a limited number of modes:

$$\mathbf{u}(t) \approx \boldsymbol{\Phi}_r \mathbf{q}(t) = \sum_{i=1}^r \boldsymbol{\phi}_i q_i(t)$$

where:
- $\boldsymbol{\Phi}_r = [\boldsymbol{\phi}_1, \boldsymbol{\phi}_2, \ldots, \boldsymbol{\phi}_r]$ = truncated mode shape matrix
- $\mathbf{q}(t) = [q_1(t), q_2(t), \ldots, q_r(t)]^T$ = modal coordinates
- $r \ll n$ = number of retained modes

### Modal Equations

Substituting the modal expansion into the governing equation and pre-multiplying by $\boldsymbol{\Phi}_r^T$:

$$\ddot{\mathbf{q}} + 2\boldsymbol{\zeta}\boldsymbol{\Omega}\dot{\mathbf{q}} + \boldsymbol{\Omega}^2\mathbf{q} = \boldsymbol{\Phi}_r^T\mathbf{f}(t)$$

where:
- $\boldsymbol{\Omega} = \text{diag}(\omega_1, \omega_2, \ldots, \omega_r)$ = natural frequency matrix
- $\boldsymbol{\zeta} = \text{diag}(\zeta_1, \zeta_2, \ldots, \zeta_r)$ = modal damping ratios

### Advantages

Modal reduction provides:
1. **Computational efficiency**: $r \ll n$ reduces system size dramatically
2. **Physical insight**: Each mode has clear interpretation
3. **Numerical stability**: Well-conditioned equations
4. **Parallelization**: Modal equations are decoupled

---

## Dynamic Simulation

### Modal Integration

Each modal equation is integrated independently:

$$\ddot{q}_i + 2\zeta_i\omega_i\dot{q}_i + \omega_i^2 q_i = \hat{f}_i(t)$$

where $\hat{f}_i(t) = \boldsymbol{\phi}_i^T\mathbf{f}(t)$ is the modal force.

### Time Integration Schemes

Several integration methods are available:

#### 1. Explicit Runge-Kutta (RK4)
4th-order accurate, conditionally stable:

```julia
# RK4 step
k1 = f(t, u)
k2 = f(t + Δt/2, u + Δt*k1/2)  
k3 = f(t + Δt/2, u + Δt*k2/2)
k4 = f(t + Δt, u + Δt*k3)
u_new = u + Δt/6 * (k1 + 2*k2 + 2*k3 + k4)
```

#### 2. Crank-Nicolson (CN)
2nd-order accurate, unconditionally stable:

$$\mathbf{u}_{n+1} = \mathbf{u}_n + \frac{\Delta t}{2}[\mathbf{f}(\mathbf{u}_n, t_n) + \mathbf{f}(\mathbf{u}_{n+1}, t_{n+1})]$$

Requires solving nonlinear system at each time step.

#### 3. Newmark-β Method
Commonly used in structural dynamics:

$$\mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t \dot{\mathbf{u}}_n + \frac{\Delta t^2}{2}[(1-2\beta)\ddot{\mathbf{u}}_n + 2\beta\ddot{\mathbf{u}}_{n+1}]$$

### Loading Patterns

#### Frequency Sweep (Chirp)
Linear frequency sweep from $f_1$ to $f_2$:

$$f(t) = F_0 \sin\left(2\pi\left[f_1 t + \frac{f_2-f_1}{2T}t^2\right]\right)$$

#### Multi-sine Excitation
Simultaneous excitation at multiple frequencies:

$$f(t) = \sum_{i=1}^n F_i \sin(2\pi f_i t + \phi_i)$$

#### Moving Load
Load progressing along bridge:

$$x_{load}(t) = v \cdot t \mod L_{bridge}$$

---

## Numerical Considerations

### Convergence Criteria

1. **Spatial convergence**: Refine mesh until modal frequencies converge
2. **Temporal convergence**: Reduce time step until solution is stable
3. **Modal convergence**: Include enough modes to capture response

### Stability Limits

- **Explicit methods**: $\Delta t < \frac{2}{\omega_{max}}$ (Courant condition)
- **Implicit methods**: Unconditionally stable but accuracy limited

### Error Sources

1. **Discretization error**: Finite element approximation
2. **Modal truncation error**: Neglected high-frequency modes  
3. **Integration error**: Time stepping approximation
4. **Round-off error**: Floating point precision

---

## References

1. **Finite Elements**: Cook, R.D. et al. "Concepts and Applications of Finite Element Analysis" (2002)
2. **Structural Dynamics**: Chopra, A.K. "Dynamics of Structures" (2017)  
3. **Modal Analysis**: Ewins, D.J. "Modal Testing: Theory, Practice and Application" (2000)
4. **Temperature Effects**: Huang, H. et al. "Temperature-dependent modal analysis of bridges" (2018)

---

**Previous**: [Getting Started](getting-started.html) | **Next**: [API Reference](api-reference.html)
