using CUDA
using LinearAlgebra

"""
    extract_modes(A, C, dt)

Extract modal parameters from state-space system identification results.

This function converts discrete-time state-space matrices from Stochastic Subspace 
Identification (SSI) to modal parameters including natural frequencies, damping ratios,
and mode shapes.

# Arguments
- `A::Matrix`: Discrete-time state transition matrix from SSI
- `C::Matrix`: Discrete-time output matrix from SSI
- `dt::Float64`: Sampling time interval used in identification

# Returns
- `freq::Vector{Float64}`: Natural frequencies in Hz
- `zeta::Vector{Float64}`: Damping ratios (dimensionless)
- `shapes_norm::Matrix{Float64}`: Normalized mode shapes

# Algorithm
1. Compute eigenvalues and eigenvectors of A matrix
2. Convert discrete-time eigenvalues to continuous-time: s = ln(λ)/dt
3. Extract natural frequencies: f = |Im(s)|/(2π)
4. Extract damping ratios: ζ = -Re(s)/|s|
5. Compute mode shapes: Ψ = C × eigenvectors
6. Normalize mode shapes to maximum absolute value

# Theory
For a discrete-time state-space system with sampling time dt:
- λ = exp(s·dt) where s are continuous-time eigenvalues
- Natural frequency: ωₙ = |Im(s)|, f = ωₙ/(2π)
- Damping ratio: ζ = -Re(s)/|s|

# Applications
- Operational Modal Analysis (OMA) from ambient vibration data
- System identification validation against FEM predictions
- Structural health monitoring parameter extraction

# See Also
- [`run_ssi_and_save`](@ref): Batch SSI processing over model orders
- [`cov_ssi_gpu`](@ref): GPU-accelerated covariance-driven SSI
"""
function extract_modes(A, C, dt)
    # Compute eigenvalues and right eigenvectors
    vals, vecs = eigen(A)
    λ = vals

    # Convert discrete-time eigenvalues to continuous-time
    s = log.(λ) ./ dt

    # Natural frequencies (Hz)
    freq = abs.(imag.(s)) ./ (2π)

    # Damping ratios
    zeta = -real.(s) ./ abs.(s)

    # Mode shapes (usually first l rows of eigenvectors, l = #outputs)
    # Normalize mode shapes (e.g. to unit norm)
    shapes = C * vecs
    shapes_norm = shapes ./ maximum(abs, shapes; dims=1)

    return freq, zeta, shapes_norm
end

"""
    run_ssi_and_save(y, i, n, dt)

Perform Stochastic Subspace Identification over multiple model orders.

Executes covariance-driven SSI for each specified model order and extracts 
modal parameters, enabling stabilization diagram analysis for optimal order selection.

# Arguments
- `y::Matrix`: Output data matrix [n_outputs × n_samples]
- `i::Int`: Number of block rows in Hankel matrices
- `n::Vector{Int}`: Vector of model orders to test
- `dt::Float64`: Sampling time interval

# Returns
- `all_freq::Vector{Vector}`: Natural frequencies for each model order
- `all_zeta::Vector{Vector}`: Damping ratios for each model order  
- `all_shapes::Vector{Matrix}`: Mode shapes for each model order

# Workflow
1. For each model order in `n`:
   - Perform covariance-driven SSI using GPU acceleration
   - Extract modal parameters using eigenvalue decomposition
   - Store results for stabilization analysis

# Usage Example
```julia
# Process ambient vibration data
model_orders = [10, 20, 30, 40, 50]
frequencies, dampings, shapes = run_ssi_and_save(response_data, 20, model_orders, 0.01)

# Analyze stabilization
for (i, order) in enumerate(model_orders)
    println("Order $order: $(length(frequencies[i])) modes identified")
end
```

# Notes
- Results should be analyzed with stabilization diagrams
- Higher model orders may identify noise modes
- GPU acceleration requires CUDA.jl and compatible hardware

# See Also
- [`extract_modes`](@ref): Modal parameter extraction from state-space matrices
- [`cov_ssi_gpu`](@ref): Core GPU SSI implementation
"""
function run_ssi_and_save(y, i, n, dt)
    all_freq = []
    all_zeta = []
    all_shapes = []

    for n_ in n
        A, C = cov_ssi_gpu(y, i, n_)
        freq, zeta, shapes = extract_modes(A, C, dt)
        push!(all_freq, freq)
        push!(all_zeta, zeta)
        push!(all_shapes, shapes)
    end

    return all_freq, all_zeta, all_shapes
end