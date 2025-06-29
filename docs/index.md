---
layout: default
title: BridgeFEM.jl - Temperature-Dependent Bridge Analysis
description: Advanced finite element analysis for bridges with temperature-dependent materials and dynamic support structures
toc: false
---

# 🌉 BridgeFEM.jl

<div class="feature-box">
<p><strong>Advanced finite element analysis for bridges with temperature-dependent materials and dynamic support structures.</strong></p>
</div>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://your-username.github.io/BridgeFEM.jl/)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()

---

## ✨ Key Features

<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1rem; margin: 2rem 0;">

<div class="feature-box">
<h3>🌡️ Temperature-Dependent Materials</h3>
<p>Realistic thermal effects with interpolated Young's modulus across temperature ranges. Perfect for seasonal analysis and thermal stress studies.</p>
</div>

<div class="feature-box">
<h3>🏗️ Dynamic Support Modeling</h3>
<p>Piers, cables, and auxiliary structures with arbitrary orientations. Fully integrated thermal expansion and support system dynamics.</p>
</div>

<div class="feature-box">
<h3>📊 Advanced Modal Analysis</h3>
<p>Eigenvalue decomposition across temperature ranges with automatic mode tracking using Modal Assurance Criterion (MAC).</p>
</div>

<div class="feature-box">
<h3>🎬 Dynamic Simulation</h3>
<p>Time-domain response with various loading patterns: chirp sweeps, multi-sine, moving loads, and custom force functions.</p>
</div>

<div class="feature-box">
<h3>📈 Rich Visualization</h3>
<p>Static plots and dynamic animations with customizable scaling. Export-ready figures and GIF animations for presentations.</p>
</div>

<div class="feature-box">
<h3>💾 Data Persistence</h3>
<p>Complete simulation configurations saved to JSON with full reproducibility. Share and version your analysis setups.</p>
</div>

</div>

---

## 🚀 Quick Example

Get started with a complete bridge analysis in just a few lines:

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
bridge = BridgeOptions(50, bc, 300.0, 7800.0, 4.0, 3.0, E_data, 50.0)

# Add vertical pier at center with same temperature dependence
pier = SupportElement(26, [1,2,3], -90.0, 5, 0.5, 0.02, E_data, 50.0, [1,2,3])

# Perform modal analysis across temperature range
temperatures = [-10.0, 20.0, 50.0]
M, K, frequencies, mode_shapes = assemble_and_decompose(bridge, temperatures, supports=[pier])

# Visualize first mode at 20°C
plot_bridge_with_supports(bridge, [pier], 
                         mode_shape=mode_shapes[:, 1, 2], 
                         scale_factor=100.0)
```

### 📊 Results Preview

The analysis reveals how bridge dynamics change with temperature:

$$f_1(-10°C) = 2.45 \text{ Hz} \quad \rightarrow \quad f_1(50°C) = 1.89 \text{ Hz}$$

**Temperature effect:** `-23%` frequency reduction from winter to summer conditions.

---

## 📚 Documentation

<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1rem; margin: 2rem 0;">

<div class="feature-box">
<h3>📖 Getting Started</h3>
<p>Installation, basic usage, and your first simulation walkthrough.</p>
<a href="getting-started.html" class="btn">Get Started →</a>
</div>

<div class="feature-box">
<h3>🔬 Theory</h3>
<p>Mathematical foundations, FEM formulation, and model order reduction techniques.</p>
<a href="theory.html" class="btn">Learn Theory →</a>
</div>

<div class="feature-box">
<h3>📋 API Reference</h3>
<p>Complete function documentation with parameters, examples, and usage notes.</p>
<a href="api-reference.html" class="btn">Browse API →</a>
</div>

<div class="feature-box">
<h3>💡 Examples</h3>
<p>Worked examples, tutorials, and case studies for real-world applications.</p>
<a href="examples.html" class="btn">See Examples →</a>
</div>

</div>

---

## 🛠️ Installation

```julia
# From Julia REPL
using Pkg
Pkg.add(url="https://github.com/your-username/BridgeFEM.jl")

# For development
Pkg.develop(url="https://github.com/your-username/BridgeFEM.jl")
```

---

## 🏗️ Core Concepts

### 🌉 Bridge Model
The `BridgeOptions` struct defines the main bridge structure with geometry, material properties, and temperature dependence. Each bridge is discretized using **2D frame elements** with 3 DOFs per node: $(u, v, \theta)$.

### 🏛️ Support Elements  
`SupportElement` represents auxiliary structures like piers or cables that can be oriented at **arbitrary angles** and have their own material properties and boundary conditions.

### 🌡️ Temperature Effects
Both bridge and supports can have **temperature-dependent stiffness**, enabling realistic thermal analysis of structures under varying environmental conditions:

$$E(T) = \text{interpolate}(T; \{T_i, E_i\}_{i=1}^n)$$

### 🎵 Modal Analysis
The library performs **eigenvalue decomposition** across temperature ranges, automatically filtering modes by frequency and normalizing mode shapes for stability.

---

## 🔬 Scientific Applications

<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1rem; margin: 2rem 0;">

<div class="feature-box">
<h3>🏗️ Structural Engineering</h3>
<ul>
<li>Bridge design verification under thermal loading</li>
<li>Support system optimization and failure analysis</li>
<li>Dynamic response prediction for various excitations</li>
</ul>
</div>

<div class="feature-box">
<h3>🔬 Research Applications</h3>
<ul>
<li>Temperature effects on structural dynamics</li>
<li>Model order reduction validation</li>
<li>System identification from ambient vibration data</li>
</ul>
</div>

<div class="feature-box">
<h3>📊 Operational Modal Analysis</h3>
<ul>
<li>GPU-accelerated Stochastic Subspace Identification</li>
<li>Automated modal parameter extraction</li>
<li>Structural health monitoring applications</li>
</ul>
</div>

</div>

---

## 📈 Architecture Overview

```
BridgeFEM.jl
├── 🏗️ Core Types
│   ├── BridgeOptions      # Main bridge definition
│   ├── SupportElement     # Auxiliary structures  
│   ├── BridgeBC          # Boundary conditions
│   └── SimulationOptions # Complete configuration
├── 🔧 Assembly
│   ├── Frame elements     # 2D beam finite elements
│   ├── Temperature deps   # Material property interpolation
│   ├── DOF mapping       # Bridge-support connectivity
│   └── Boundary conds    # Constraint enforcement
├── 🎵 Modal Analysis
│   ├── Eigenvalue solver  # ARPACK integration
│   ├── Mode tracking     # Modal Assurance Criterion
│   └── Normalization     # Mass and length scaling
├── 🎬 Dynamic Simulation
│   ├── Modal reduction    # Frequency-based truncation
│   ├── Time integration   # DifferentialEquations.jl
│   └── Temperature interp # Real-time property updates
└── 📊 Visualization
    ├── Static plots       # Mode shapes and geometry
    ├── Animations         # Dynamic response movies
    └── Export formats     # PNG, GIF, PDF output
```

---

## 🤝 Contributing

Contributions are welcome! Please see our [contribution guidelines](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone repository
git clone https://github.com/your-username/BridgeFEM.jl.git
cd BridgeFEM.jl

# Set up Julia environment
julia --project=. -e "using Pkg; Pkg.instantiate()"

# Run tests
julia --project=. -e "using Pkg; Pkg.test()"
```

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

<div class="text-center" style="margin: 3rem 0;">
<p><strong>Ready to analyze your bridge?</strong></p>
<a href="getting-started.html" class="btn" style="font-size: 1.2em; padding: 1rem 2rem;">Get Started Now →</a>
</div>
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
