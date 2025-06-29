# BridgeFEM.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://jgrashorn.github.io/BridgeFEM.jl/)

A Julia library for finite element analysis of bridges with temperature-dependent materials, dynamic support structures, and advanced modal analysis capabilities.

## ✨ Features

- **🌡️ Temperature-dependent materials**: Realistic thermal effects with interpolated Young's modulus
- **🏗️ Dynamic support modeling**: Piers, cables, and auxiliary structures with arbitrary orientations  
- **📊 Advanced modal analysis**: Eigenvalue decomposition across temperature ranges with mode tracking
- **🎬 Dynamic simulation**: Time-domain response with various loading patterns
- **📈 Comprehensive visualization**: Static plots and dynamic animations with deformation scaling
- **💾 Data persistence**: Complete simulation configurations saved to JSON
- **🔍 Operational Modal Analysis**: GPU-accelerated system identification from ambient data

## Disclaimer!

As is evident by the heavy use of emojis in the readme and throughout the examples, many parts of this code (and nearly all of the documentation) were created using LLMs. The core functionality should work as intended, but there might be errors in the explanations.

## 🚀 Quick Start

### Installation

```julia
# From Julia REPL
pkg> add https://github.com/jgrashorn/BridgeFEM.jl
```

### Basic Example

```julia
using BridgeFEM

# Define temperature-dependent steel properties
E_data = [
    -10.0  250e9;   # Young's modulus at -10°C (winter)
     20.0  207e9;   # Young's modulus at 20°C (reference)  
     50.0  150e9    # Young's modulus at 50°C (summer)
]

# Create boundary conditions: fixed-simply supported
bc = BridgeBC([[1, "all"], [51, "y"]])

# Define 300m bridge with 50 elements
bridge = BridgeOptions(
    50,           # Number of elements
    bc,           # Boundary conditions
    300.0,        # Length (m)
    7800.0,       # Density (kg/m³)
    4.0,          # Cross-sectional area (m²)
    3.0,          # Moment of inertia (m⁴)
    E_data,       # Temperature-E data
    50.0          # Cutoff frequency (Hz)
)

# Add vertical pier at center
pier = SupportElement(
    26,           # Connection node
    [1,2,3],      # Connect all DOFs
    -90.0,        # Vertical downward
    5,            # 5 elements
    0.5,          # Area (m²)
    0.02,         # Inertia (m⁴)
    E_data,       # Same material as bridge
    50.0,         # Height (m)
    [1,2,3]       # Fixed base
)

# Perform modal analysis
temperatures = [-10.0, 20.0, 50.0]
M, K, frequencies, mode_shapes = assemble_and_decompose(
    bridge, temperatures, supports=[pier]
)

# Visualize first mode at 20°C
plot_bridge_with_supports(bridge, [pier], 
                         mode_shape=mode_shapes[:, 1, 2], 
                         scale_factor=100.0)
```

### Dynamic Analysis

```julia
# Create frequency sweep loading
load_func = create_chirp_loading(0.1, 10.0, target_node=26)

# Solve dynamic response
sol = solve_dynamics(bridge, [pier], temperatures, load_func)

# Create animation
animate_dynamic_response(bridge, [pier], sol, scale_factor=50.0)
```

## 📚 Documentation

Complete documentation is available at: **[BridgeFEM.jl Documentation](https://jgrashorn.github.io/BridgeFEM.jl/)**

### Documentation Sections

- **[Getting Started](docs/getting-started.md)**: Installation and first simulation
- **[Theory](docs/theory.md)**: FEM formulation, model order reduction, and mathematical background
- **[API Reference](docs/api-reference.md)**: Complete function and type documentation
- **[Examples](examples/)**: Worked examples and case studies

## Scientific Applications

### Structural Engineering
- Bridge design verification under thermal loading
- Dynamic response prediction for various excitation types
- Support system optimization and failure analysis

### Research Applications  
- Temperature effects on structural dynamics
- Model order reduction validation
- System identification from ambient vibration data

### Operational Modal Analysis
- GPU-accelerated Stochastic Subspace Identification (SSI)
- Automated modal parameter extraction
- Structural health monitoring applications

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```julia
# Clone repository
git clone https://github.com/jgrashorn/BridgeFEM.jl.git
cd BridgeFEM.jl

# Activate and install dependencies
julia --project=. -e "using Pkg; Pkg.instantiate()"

# Run tests
julia --project=. -e "using Pkg; Pkg.test()"
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use BridgeFEM.jl in your research, please cite:

```bibtex
@software{bridgefem_jl,
    title = {BridgeFEM.jl: Temperature-dependent Bridge Analysis with Julia},
    author = {Jan Grashorn},
    year = {2025},
    url = {https://github.com/jgrashorn/BridgeFEM.jl}
}
```
