# Quimb

## Official Resources
- **Homepage**: https://quimb.readthedocs.io/
- **Repository**: https://github.com/jcmgray/quimb
- **License**: Apache License 2.0
- **PyPI**: https://pypi.org/project/quimb/

## Overview
Quimb is a standalone Python library for quantum information and many-body calculations. It aims to bridge the gap between "easy" and "fast" by providing a high-level, interactive API for defining quantum states and operators, coupled with a powerful tensor network engine (`quimb.tensor`). Quimb handles arbitrary tensor network geometries (1D, 2D, 3D, hyperbolic) and integrates with advanced contraction path optimizers and hardware-accelerated backends (like CuPy, JAX, and Torch) to perform simulations of large-scale quantum systems and quantum circuits efficiently.

**Scientific domain**: Quantum Information, Many-Body Physics, Quantum Circuit Simulation.
**Target user community**: Researchers in quantum computing and condensed matter physics needing rapid prototyping and flexible geometry support.

## Theoretical Methods
- **Tensor Network Algorithms**: DMRG, TEBD, simple-update PEPS, HOTRG.
- **Quantum Information**: Entanglement entropy, negativity, fidelity, partial trace.
- **Dynamics**: Real-time evolution using TEBD, TDVP (via extensions).
- **Quantum Circuits**: Circuit simulation via tensor network contraction.

## Capabilities (CRITICAL)
- **Geometry Agnostic**: Natively handles MPS, PEPS, MERA, and random tensor graphs.
- **Advanced Contraction**: Uses `cotengra` (COntraction Tree ENGRAne) to find optimal contraction paths for complex networks.
- **Backend Acceleration**: Can dispatch tensor operations to GPU-enabled backends (JAX, TensorFlow, PyTorch).
- **Visualization**: Excellent drawing tools for tensor networks, showing bond dimensions and connectivity.
- **Interactive**: Designed for Jupyter notebooks with rich visual outputs.
- **Exact Diagonalization**: Includes tools for full Hilbert space caching and ED.

## Key Strengths
### Flexibility
- Not restricted to lattices; can handle random graphs or quantum circuits with arbitrary connectivity.

### Performance
- "Easy but Fast": Uses Numba for critical sections.
- Supports distributed computing via `ray` or `dask`.

### Visualization
- Best-in-class visualization for tensor networks, greatly aiding debugging and intuition.

## Inputs & Outputs
- **Input**: Quantum circuits (qasm), Hamiltonians, Tensor definitions.
- **Output**: Entanglement measures, Expectation values, Optimized states, Visualizations (.png, .svg).

## Interfaces & Ecosystem
- **Python**: Pure Python installation.
- **Integration**: Works with `cotengra` for path finding, `slepc4py` for eigensolvers.
- **Backends**: Compatible with `numpy`, `jax`, `torch`, `tensorflow`, `cupy`.

## Advanced Features
- **Quantum Circuit Simulation**: Simulates deep circuits by optimizing contraction order.
- **Hyper-optimization**: Can optimize contraction trees for repeated execution.
- **Symmetries**: Basic support for symmetries.

## Performance Characteristics
- **Contraction**: State-of-the-art contraction path finding makes it competitive for "hard" networks (e.g., Sycamore circuits).
- **Speed**: GPU backends provide massive speedups for large bond dimensions.

## Comparisons with Other Codes
- **vs TeNPy**: TeNPy is more specialized for 1D/2D ground states with U(1) symmetry and physics-heavy features. Quimb is better for random circuits, visualization, and geometric flexibility.
- **vs TensorNetwork**: Quimb is a higher-level framework that *can* use TensorNetwork as a backend. Quimb provides the physics/QI logic, TN provides the contraction engine.
- **vs Qiskit**: Quimb simulates circuits using tensor networks, often scaling better for low-entanglement circuits than state-vector simulators.

## Application Areas
- **Quantum Chaos**: Studying entanglement growth in random circuits.
- **NISQ Simulation**: Simulating Google/IBM quantum supremacy circuits.
- **Tensor Network Machine Learning**: Training generative models.

## Best Practices
- **Path Optimization**: Always pre-optimize contraction paths for large networks (reuse the path).
- **Backends**: Switch to `jax` or `cupy` for GPU execution.
- **Visualization**: Use `.draw()` frequently to inspect network structure.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/jcmgray/quimb
2. Documentation: https://quimb.readthedocs.io/
3. Papers citing Quimb (e.g., Gray & Kourtis, arXiv:2002.01935).

**Confidence**: VERIFIED - Active development and widely used.
**Verification status**: ✅ VERIFIED
