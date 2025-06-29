---
layout: default
title: Theory - BridgeFEM.jl
description: Mathematical foundations of temperature-dependent finite element analysis
toc: true
---

# 📐 Theoretical Background

This page describes the mathematical foundations, finite element formulation, and model order reduction techniques used in BridgeFEM.jl.

> **🎯 Key Topics**: Finite Element Method • Temperature-Dependent Materials • Modal Analysis • Model Order Reduction • Dynamic Simulation

---

## 🔬 Finite Element Method

### Governing Equations

For a structural system, the equations of motion in physical coordinates are:

$$\boxed{\mathbf{M}\ddot{\mathbf{u}} + \mathbf{C}\dot{\mathbf{u}} + \mathbf{K}(T)\mathbf{u} = \mathbf{f}(t)}$$

**System Matrices:**
- $\mathbf{M} \in \mathbb{R}^{n \times n}$ = global mass matrix
- $\mathbf{C} \in \mathbb{R}^{n \times n}$ = damping matrix  
- $\mathbf{K}(T) \in \mathbb{R}^{n \times n}$ = temperature-dependent stiffness matrix
- $\mathbf{u}(t) \in \mathbb{R}^n$ = displacement vector
- $\mathbf{f}(t) \in \mathbb{R}^n$ = external force vector

### Discretization Strategy

The bridge structure is discretized using **2D frame elements** with 3 degrees of freedom per node:

$$\mathbf{u}_i = \begin{bmatrix} u_i \\ v_i \\ \theta_i \end{bmatrix} \quad \text{where} \quad \begin{cases}
u_i & = \text{horizontal displacement} \\
v_i & = \text{vertical displacement} \\
\theta_i & = \text{rotation about z-axis}
\end{cases}$$

**System Size:** For $n$ elements with $n+1$ nodes → **Total DOFs = $3(n+1)$**

---

## 🏗️ 2D Frame Elements

### Element Stiffness Matrix

Each frame element connects two nodes with 6 total DOFs. The **local element stiffness matrix** is:

$$\mathbf{k}^e = \begin{bmatrix}
\dfrac{EA}{L} & 0 & 0 & -\dfrac{EA}{L} & 0 & 0 \\[0.3em]
0 & \dfrac{12EI}{L^3} & \dfrac{6EI}{L^2} & 0 & -\dfrac{12EI}{L^3} & \dfrac{6EI}{L^2} \\[0.3em]
0 & \dfrac{6EI}{L^2} & \dfrac{4EI}{L} & 0 & -\dfrac{6EI}{L^2} & \dfrac{2EI}{L} \\[0.3em]
-\dfrac{EA}{L} & 0 & 0 & \dfrac{EA}{L} & 0 & 0 \\[0.3em]
0 & -\dfrac{12EI}{L^3} & -\dfrac{6EI}{L^2} & 0 & \dfrac{12EI}{L^3} & -\dfrac{6EI}{L^2} \\[0.3em]
0 & \dfrac{6EI}{L^2} & \dfrac{2EI}{L} & 0 & -\dfrac{6EI}{L^2} & \dfrac{4EI}{L}
\end{bmatrix}$$

**Material and Geometric Parameters:**
- $E(T)$ = Temperature-dependent Young's modulus [Pa]
- $A$ = Cross-sectional area [m²]
- $I$ = Second moment of inertia [m⁴]
- $L$ = Element length [m]

### Element Mass Matrix

The **consistent mass matrix** derived from Hermite interpolation functions:

$$\mathbf{m}^e = \frac{\rho A L}{420} \begin{bmatrix}
140 & 0 & 0 & 70 & 0 & 0 \\
0 & 156 & 22L & 0 & 54 & -13L \\
0 & 22L & 4L^2 & 0 & 13L & -3L^2 \\
70 & 0 & 0 & 140 & 0 & 0 \\
0 & 54 & 13L & 0 & 156 & -22L \\
0 & -13L & -3L^2 & 0 & -22L & 4L^2
\end{bmatrix}$$

**Properties:**
- $\rho$ = Material density [kg/m³]
- Preserves total element mass: $\int_0^L \rho A \, dx = \rho A L$
- Couples translational and rotational inertia

### Assembly Process

**Global matrix assembly** using standard FEM connectivity:

```julia
function assemble_matrices(bridge, T)
    M = spzeros(n_dofs, n_dofs)
    K = spzeros(n_dofs, n_dofs)
    
    for e = 1:n_elements
        # Temperature-dependent element matrices
        ke = frame_elem_stiffness(E(T)*A, E(T)*I, L_element)
        me = frame_elem_mass(ρ, A, L_element)
        
        # DOF connectivity mapping
        dofs = [3*(e-1)+1:3*(e-1)+3; 3*e+1:3*e+3]
        
        # Assemble into global matrices
        K[dofs, dofs] += ke
        M[dofs, dofs] += me
    end
    
    return M, K
end
```

---

## 🌡️ Temperature Dependence

### Material Property Interpolation

Young's modulus varies with temperature according to **user-defined data points**:

$$\boxed{E(T) = \text{interpolate}(T; \{(T_i, E_i)\}_{i=1}^{n_{points}})}$$

**Implementation using linear interpolation:**

```julia
# Define temperature-Young's modulus pairs
E_data = [
    -10.0  250e9;   # Winter conditions
     20.0  207e9;   # Reference temperature  
     50.0  150e9    # Summer conditions
]

# Create interpolation function
E_interp = interpolate((E_data[:,1],), E_data[:,2], Gridded(Linear()))
E = T -> E_interp(T)  # Temperature-dependent function
```

### Stiffness Matrix Evolution

At temperature $T$, the **global stiffness matrix** becomes:

$$\mathbf{K}(T) = \sum_{e=1}^{n_{elem}} \mathbf{L}_e^T \mathbf{k}^e(E(T)) \mathbf{L}_e$$

where:
- $\mathbf{L}_e \in \{0,1\}^{n \times 6}$ = Boolean connectivity matrix for element $e$
- $\mathbf{k}^e(E(T))$ = Element stiffness matrix at temperature $T$

**Key Properties:**
- Mass matrix $\mathbf{M}$ remains **temperature-independent**
- Only stiffness scales with $E(T)$: $\mathbf{k}^e(T) = E(T) \cdot \mathbf{k}^e_0$

---

## 🏛️ Support Modeling

### Coordinate Transformation

Support elements (piers, cables, etc.) can be oriented at **arbitrary angles**. The transformation from local to global coordinates is:

$$\boxed{\mathbf{k}_{global}^{support} = \mathbf{T}^T \mathbf{k}_{local}^{support} \mathbf{T}}$$

The **2D rotation transformation matrix** for angle $\theta$ (measured counterclockwise from horizontal):

$$\mathbf{T} = \begin{bmatrix}
\cos\theta & \sin\theta & 0 & 0 & 0 & 0 \\
-\sin\theta & \cos\theta & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & \cos\theta & \sin\theta & 0 \\
0 & 0 & 0 & -\sin\theta & \cos\theta & 0 \\
0 & 0 & 0 & 0 & 0 & 1
\end{bmatrix}$$

**Angle Conventions:**
- $\theta = 0°$ → Horizontal, pointing right
- $\theta = -90°$ → Vertical, pointing down  
- $\theta = 90°$ → Vertical, pointing up

### DOF Connectivity

Support elements connect to the main bridge through **shared degrees of freedom**:

1. **Bridge DOFs**: $\{u_i, v_i, \theta_i\}$ at connection node $i$
2. **Support DOFs**: Local support numbering → mapped to global system  
3. **Connectivity Types**:
   - **Pinned**: $[1,2]$ → Share only translations
   - **Rigid**: $[1,2,3]$ → Share all DOFs
   - **Custom**: User-defined DOF sharing

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

**System Properties:**
- **Expanded DOFs**: $n_{total} = n_{bridge} + \sum n_{support,i}$  
- **Coupling**: Through shared DOF constraints
- **Temperature dependence**: Both bridge and supports scale with $E(T)$

---

## 🎵 Modal Analysis

### Generalized Eigenvalue Problem

Natural frequencies and mode shapes are obtained by solving:

$$\boxed{\mathbf{K}(T)\boldsymbol{\phi}_i = \lambda_i(T) \mathbf{M}\boldsymbol{\phi}_i}$$

**Solution Properties:**
- $\lambda_i(T) = \omega_i^2(T)$ = squared natural frequency
- $\boldsymbol{\phi}_i(T) \in \mathbb{R}^n$ = mode shape vector
- $f_i(T) = \frac{\omega_i(T)}{2\pi}$ = natural frequency [Hz]

### Temperature-Dependent Eigenvalues

Since $\mathbf{K}(T) = E(T) \mathbf{K}_0$, the **frequencies scale as**:

$$f_i(T) = \sqrt{\frac{E(T)}{E_0}} \cdot f_{i,0}$$

This relationship enables **efficient frequency interpolation** across temperature ranges.

### Solution Method

For large sparse systems, **ARPACK** provides efficient computation:

```julia
# Solve for lowest n_modes eigenvalues
λ, φ = eigs(K, M, nev=n_modes, which=:SM, sigma=1e-6)

# Convert to natural frequencies  
frequencies = sqrt.(real.(λ)) / (2π)  # [Hz]
mode_shapes = real.(φ)               # [n_dofs × n_modes]
```

### Mass Normalization

Mode shapes are normalized to **unit modal mass**:

$$\boldsymbol{\phi}_i^T \mathbf{M} \boldsymbol{\phi}_i = 1 \quad \forall i$$

This ensures the modal equations have **consistent scaling**.

---

## 📉 Model Order Reduction

### Modal Truncation Strategy

Physical displacements are approximated using a **reduced basis**:

$$\boxed{\mathbf{u}(t) \approx \boldsymbol{\Phi}_r \mathbf{q}(t) = \sum_{i=1}^r \boldsymbol{\phi}_i q_i(t)}$$

**Truncation Criteria:**
- **Frequency-based**: Retain modes up to cutoff frequency $f_{cut}$
- **Number-based**: Keep first $r$ modes
- **Energy-based**: Retain modes capturing desired energy fraction

### Modal Coordinate Equations

Substituting modal expansion and **pre-multiplying by $\boldsymbol{\Phi}_r^T$**:

$$\boxed{\ddot{\mathbf{q}} + 2\boldsymbol{\zeta}\boldsymbol{\Omega}\dot{\mathbf{q}} + \boldsymbol{\Omega}^2\mathbf{q} = \boldsymbol{\Phi}_r^T\mathbf{f}(t)}$$

**Modal System Matrices:**
- $\boldsymbol{\Omega}(T) = \text{diag}(\omega_1(T), \ldots, \omega_r(T))$ = natural frequency matrix
- $\boldsymbol{\zeta} = \text{diag}(\zeta_1, \ldots, \zeta_r)$ = modal damping ratios
- $\mathbf{q}(t) \in \mathbb{R}^r$ = modal coordinates

### Computational Advantages

Modal reduction provides:

1. **Dimension reduction**: $n \to r$ where $r \ll n$
2. **Decoupled equations**: Diagonal system matrices
3. **Physical insight**: Each mode represents distinct vibration pattern
4. **Computational efficiency**: $O(r^3)$ vs $O(n^3)$ for time integration

---

## 🎬 Dynamic Simulation

### Temperature-Dependent Modal Parameters

For time-varying temperature $T(t)$, modal parameters are **interpolated**:

$$\omega_i(t) = \omega_i(T(t)), \quad \boldsymbol{\phi}_i(t) = \boldsymbol{\phi}_i(T(t))$$

**Implementation using interpolation functions:**

```julia
# Create interpolators for each parameter
ω_interp = interpolate((modes, temps), frequencies, Gridded(Linear()))
Φ_interp = interpolate((dofs, modes, temps), mode_shapes, Gridded(Linear()))

# Time-dependent functions
ω_T = t -> ω_interp(1:n_modes, T(t))
Φ_T = t -> Φ_interp(1:n_dofs, 1:n_modes, T(t))
```

### Modal ODE System

The **temperature-dependent modal equations** become:

$$\ddot{q}_i + 2\zeta_i\omega_i(T(t))\dot{q}_i + \omega_i^2(T(t))q_i = \hat{f}_i(t)$$

where $\hat{f}_i(t) = \boldsymbol{\phi}_i^T(T(t)) \mathbf{f}(t)$ is the **modal force**.

### Physical Reconstruction

Modal coordinates are transformed back to **physical displacements**:

$$\mathbf{u}(t) = \boldsymbol{\Phi}(T(t)) \mathbf{q}(t)$$

**Velocity and acceleration** are obtained by differentiation:
$$\dot{\mathbf{u}}(t) = \boldsymbol{\Phi}(T(t)) \dot{\mathbf{q}}(t) + \dot{\boldsymbol{\Phi}}(T(t)) \mathbf{q}(t)$$

---

## 📊 Modal Assurance Criterion

### Mode Tracking Across Temperature

To ensure **consistent mode ordering** across temperatures, the Modal Assurance Criterion (MAC) is used:

$$\text{MAC}_{ij} = \frac{|\boldsymbol{\phi}_i^T \boldsymbol{\phi}_j|^2}{(\boldsymbol{\phi}_i^T \boldsymbol{\phi}_i)(\boldsymbol{\phi}_j^T \boldsymbol{\phi}_j)}$$

**Interpretation:**
- $\text{MAC} = 1$ → Perfect correlation (same mode)
- $\text{MAC} = 0$ → Orthogonal modes (different modes)
- $\text{MAC} > 0.9$ → Strong correlation (likely same mode)

### Implementation

```julia
function track_modes(φ_prev, φ_curr)
    for i = 1:n_modes
        # Check correlation with previous temperature
        correlation = dot(φ_curr[:, i], φ_prev[:, i])
        
        # Flip mode if negatively correlated
        if correlation < 0
            φ_curr[:, i] *= -1
        end
    end
    return φ_curr
end
```

---

## 📚 References

1. **Bathe, K.J.** (2014). *Finite Element Procedures*. Prentice Hall.
2. **Chopra, A.K.** (2017). *Dynamics of Structures: Theory and Applications*. Pearson.
3. **Zienkiewicz, O.C. & Taylor, R.L.** (2000). *The Finite Element Method*. Butterworth-Heinemann.
4. **Allemang, R.J.** (2003). "The Modal Assurance Criterion - Twenty years of use and abuse." *Sound and Vibration*, 37(8), 14-21.
5. **Peeters, B. & De Roeck, G.** (1999). "Reference-based stochastic subspace identification for output-only modal analysis." *Mechanical Systems and Signal Processing*, 13(6), 855-878.

---

**🔬 Next Steps:** [Getting Started Guide](getting-started.md) • [API Reference](api-reference.md) • [Examples](examples.md)
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
