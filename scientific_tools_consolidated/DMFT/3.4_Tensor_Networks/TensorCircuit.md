# TensorCircuit

## Official Resources
- **Homepage**: https://tensorcircuit.readthedocs.io/
- **Repository**: https://github.com/tencent-quantum-lab/tensorcircuit
- **License**: Apache License 2.0
- **Developer**: Tencent Quantum Lab

## Overview
TensorCircuit is a next-generation open-source quantum software framework developed by Tencent. It is built on a high-performance tensor network simulation engine and is fully compatible with modern Deep Learning (DL) frameworks like TensorFlow, JAX, and PyTorch. TensorCircuit is designed for the Noisy Intermediate-Scale Quantum (NISQ) era, enabling efficient simulation of large-scale quantum circuits, variational quantum algorithms (VQA), and quantum machine learning (QML) tasks by leveraging automatic differentiation, JIT compilation, and vectorization.

**Scientific domain**: Quantum Computing, Quantum Machine Learning, Tensor Networks.
**Target user community**: Quantum algorithm researchers, QML practitioners.

## Theoretical Methods
- **Tensor Network Contraction**: Simulates quantum circuits by contracting the corresponding tensor network.
- **Automatic Differentiation (AD)**: Backpropagation through quantum circuits for gradient-based optimization.
- **Vectorization**: Parallel evaluation of batched circuits (VMAP).
- **JIT Compilation**: Just-In-Time compilation for high-speed execution.

## Capabilities (CRITICAL)
- **Circuit Simulation**: Simulates noise-free and noisy quantum circuits (density matrix, Monte Carlo trajectories).
- **Backends**: Seamless switching between TensorFlow, JAX, PyTorch, and NumPy.
- **Gradient Calculation**: Exact gradients via AD, avoiding parameter shaft limitations of hardware.
- **Visualization**: Built-in circuit and tensor network visualization.
- **Qiskit Interop**: Convert Qiskit circuits to TensorCircuit.
- **Advanced Features**: Fermion Gaussian states, specific noise models, error mitigation (ZNE, RC).

## Key Strengths
### Performance
- **Speed**: 10x-1000x faster than standard simulators (Qiskit/Cirq) for deep circuits via TN contraction engine.
- **Memory**: Tensor networks can be more memory-efficient than state vectors for certain low-entanglement circuits.

### ML Integration
- Native integration with differentiable programming makes it ideal for Variational Quantum Eigensolvers (VQE) and QNNs.

## Inputs & Outputs
- **Input**: Quantum circuits (native API or Qiskit/OpenQASM), Hamiltonians.
- **Output**: Density matrices, Expectation values, Gradients, Samples.

## Interfaces & Ecosystem
- **Ecosystem**: Works with standard Python ML stack.
- **Hardware**: Can deploy circuits to real quantum hardware (experimental).
- **Cloud**: Integration with Tencent Cloud (if applicable).

## Advanced Features
- **Parallelization**: Multi-GPU support for batched circuit evaluation.
- **Contraction Optimization**: Uses `cotengra` to find optimal contraction paths.

## Performance Characteristics
- **Benchmarks**: Outperforms TensorFlow Quantum and PennyLane in many VQA benchmarks.
- **Scalability**: Can simulate hundreds of qubits for shallow/structured circuits (e.g., 600+ qubit VQE 1D).

## Computational Cost
- **High Entanglement**: Cost grows exponentially with entanglement "treewidth", similar to other TN codes.
- **Low Entanglement**: Extremely efficient.

## Limitations & Known Constraints
- **Complexity**: Deep, highly entangled circuits (high treewidth) eventually hit the memory wall.
- **Backend Dependency**: Performance depends heavily on the chosen backend (JAX usually fastest).

## Comparison with Other Codes
- **vs Qiskit**: TensorCircuit uses TN contraction (better for some large/shallow circuits) vs Qiskit's statevector (exponential memory). TC supports native AD.
- **vs PennyLane**: Similar ML integration, but TC is built specifically on a tensor network engine, often offering better simulation performance for VQAs.
- **vs Google TensorNetwork**: TensorCircuit is a higher-level framework *built* using TN concepts (and potentially libraries), focused specifically on quantum circuits.

## Application Areas
- **VQE/QAOA**: Variational algorithms for chemistry and optimization.
- **Quantum Machine Learning**: Training quantum neural networks.
- **Error Mitigation**: Researching zero-noise extrapolation and other techniques.

## Best Practices
- **JIT**: Always use JIT (`@tc.jit`) for production runs.
- **Batching**: Use `vmap` for evaluating parameterized circuits in parallel.
- **Backend**: JAX is recommended for the best performance.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/tencent-quantum-lab/tensorcircuit
2. Paper: "TensorCircuit: A Quantum Software Framework for the NISQ Era" (arXiv:2205.10091).

**Confidence**: VERIFIED - Active maintainer (Tencent), benchmarking papers available.
**Verification status**: ✅ VERIFIED
