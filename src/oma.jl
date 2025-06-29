using CUDA
using LinearAlgebra

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