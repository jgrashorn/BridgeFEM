---
layout: default
title: BridgeFEM.jl
---

# BridgeFEM.jl

A Julia library for finite element analysis of bridges with temperature-dependent materials and dynamic, arbitrary-ish support structures (currently only straight supports are implemented).

## Quick Example

```julia
using BridgeFEM

# Create a 300m bridge with temperature-dependent steel
E_data = [
    -10.0  250e9;   # E at -10°C (winter)
     20.0  207e9;   # E at 20°C (reference)  
     50.0  150e9    # E at 50°C (summer)
]

bc = BridgeBC([[1, "all"], [51, "y"]])  # Fixed-simply supported
bridge = BridgeOptions(50, bc, 300.0, 7800.0, 4.0, 3.0, E_data, 50.0)

# Add a vertical pier at the center
pier = SupportElement(26, [1,2,3], -90.0, 5, 0.5, 0.02, E_data, 50.0, [1,2,3])

# Perform modal analysis across temperature range
temperatures = [-10.0, 20.0, 50.0]
M, K, frequencies, mode_shapes = assemble_and_decompose(bridge, temperatures, supports=[pier])

# Visualize first mode at 20°C
plot_bridge_with_supports(bridge, [pier], 
                         mode_shape=mode_shapes[:, 1, 2], 
                         scale_factor=100.0)

# Dynamic response to frequency sweep
load_vector = create_chirp_loading(0.1, 10.0, target_node=26)
sol = solve_dynamics(bridge, [pier], temperatures, load_vector)

# Create animation
animate_dynamic_response(bridge, [pier], sol.u, sol.t, 
                        scale_factor=5000.0, 
                        filename="bridge_response.gif")
```

## Getting Started

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://jgrashorn.github.io/BridgeFEM.jl/)
[![Build Status](https://github.com/jgrashorn/BridgeFEM.jl/workflows/CI/badge.svg)](https://github.com/jgrashorn/BridgeFEM.jl/actions)

### Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jgrashorn/BridgeFEM.jl")
```

### Quick Links

- **[Getting Started](getting-started.html)** - Installation, basic usage, and first simulation
- **[Theory](theory.html)** - Finite element formulation and model order reduction
- **[API Reference](api-reference.html)** - Detailed function documentation  
- **[Examples](examples.html)** - Complete worked examples and tutorials

## Key Concepts

### Bridge Model
The `BridgeOptions` struct defines the main bridge structure with geometry, material properties, and temperature dependence. Each bridge is discretized using 2D frame elements with 3 DOFs per node (u, v, θ).

### Support Elements  
`SupportElement` represents auxiliary structures like piers or cables that can be oriented at arbitrary angles and have their own material properties and boundary conditions.

### Temperature Effects
Both bridge and supports can have temperature-dependent stiffness, enabling realistic thermal analysis of structures under varying environmental conditions.

### Modal Analysis
The library performs eigenvalue decomposition across temperature ranges, automatically filtering modes by frequency and normalizing mode shapes for stability.

## Applications

- **Bridge design verification**: Modal analysis under different thermal conditions
- **Seismic analysis**: Dynamic response to ground motion with realistic support modeling  
- **Health monitoring**: Comparison of measured vs. predicted modal properties
- **Parametric studies**: Effect of support stiffness, temperature, and loading on dynamic behavior

## Architecture

```
BridgeFEM.jl
├── Core Types
│   ├── BridgeOptions      # Main bridge definition
│   ├── SupportElement     # Auxiliary structures
│   └── BridgeBC          # Boundary conditions
├── Assembly
│   ├── Frame elements     # 2D beam finite elements
│   ├── Temperature deps   # Material property interpolation
│   └── DOF mapping       # Bridge-support connectivity
├── Analysis  
│   ├── Modal analysis     # Eigenvalue problems
│   ├── Dynamic simulation # Time integration
│   └── Model reduction    # Modal truncation
└── Visualization
    ├── Static plots       # Structure and mode shapes
    └── Animations         # Dynamic response
```

## Performance

The library is optimized for efficiency:
- **Sparse matrices** for large DOF systems
- **ARPACK eigensolvers** for modal analysis
- **Interpolated properties** for temperature dependence
- **Manual integrators** available for long simulations

Typical performance on a modern laptop:
- **Assembly**: ~0.1s for 1000 DOF system
- **Modal analysis**: ~1s for first 50 modes  
- **Dynamic simulation**: ~10s for 100s at 0.01s timestep
- **Animation**: ~30s for 100 frame GIF

## Citation

If you use this library in research, please cite:

```bibtex
@software{bridgefem2025,
  title = {BridgeFEM.jl: Finite Element Analysis of Bridges with Temperature-Dependent Materials},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/jgrashorn/BridgeFEM.jl},
  doi = {10.5281/zenodo.XXXXXXX}
}
```

## License

MIT License - see [LICENSE](https://github.com/jgrashorn/BridgeFEM.jl/blob/main/LICENSE) file for details.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](https://github.com/jgrashorn/BridgeFEM.jl/blob/main/CONTRIBUTING.md) for guidelines.

---

**Next**: [Getting Started](getting-started.html) | **See also**: [Examples](examples.html)
