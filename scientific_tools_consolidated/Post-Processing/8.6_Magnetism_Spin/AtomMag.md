# AtomMag

## Official Resources
- Source Repository: https://github.com/jhu238/AtomMag
- License: Open source

## Overview
**AtomMag** is a GPU-parallel atomistic spin dynamics model developed at the University of Wisconsin-Madison. It achieves 66x speedup over CPU implementations and is validated against fidimag and analytical results, supporting large-scale atomistic magnetic simulations.

**Scientific domain**: GPU-accelerated atomistic spin dynamics  
**Target user community**: Researchers needing fast GPU atomistic spin dynamics for large magnetic systems

## Theoretical Methods
- Atomistic Landau-Lifshitz-Gilbert (LLG) equation
- Heisenberg Hamiltonian
- GPU parallelization (CUDA)
- Monte Carlo simulation
- Exchange, anisotropy, Zeeman energies
- Dzyaloshinskii-Moriya interaction

## Capabilities (CRITICAL)
- GPU-accelerated atomistic spin dynamics
- 66x speedup over CPU
- Large-scale spin systems
- Heisenberg model simulation
- Monte Carlo for equilibrium
- Validated against fidimag and analytical results

**Sources**: GitHub repository

## Key Strengths

### GPU Performance:
- 66x speedup over CPU
- CUDA parallelization
- Large system sizes feasible
- Efficient memory access

### Validation:
- Compared with fidimag
- Analytical result validation
- Reproducible benchmarks
- Well-tested

## Inputs & Outputs
- **Input formats**:
  - Configuration files
  - Exchange parameters
  - Material parameters
  
- **Output data types**:
  - Magnetization vs time
  - Energy vs time
  - Spin configurations
  - Statistical averages

## Interfaces & Ecosystem
- **CUDA**: GPU acceleration
- **C/C++**: Core computation
- **Python**: Post-processing

## Performance Characteristics
- **Speed**: Very fast (GPU)
- **Accuracy**: Good (validated)
- **System size**: Millions of spins
- **Memory**: GPU memory limited

## Computational Cost
- **Small systems**: Seconds
- **Large systems**: Minutes
- **Typical**: Very efficient

## Limitations & Known Constraints
- **NVIDIA GPU only**: Requires CUDA
- **Limited documentation**: Research code
- **Limited features**: Basic Heisenberg model
- **Small community**: Research group code

## Comparison with Other Codes
- **vs fidimag**: AtomMag is GPU (66x faster), fidimag is CPU
- **vs UppASD**: AtomMag is GPU, UppASD is CPU with more features
- **vs Spirit**: AtomMag is GPU-only, Spirit has GUI
- **Unique strength**: GPU-parallel atomistic spin dynamics with 66x speedup, validated

## Application Areas

### Large-Scale Spin Dynamics:
- Extended magnetic systems
- Long-time dynamics
- Parameter sweeps
- High-throughput simulation

### Benchmarking:
- GPU vs CPU comparison
- Code validation
- Performance testing
- Method verification

## Best Practices

### GPU Setup:
- Use NVIDIA GPU
- Monitor GPU memory
- Optimize block/grid sizes
- Compare with CPU for validation

### Simulation:
- Validate against analytical results
- Use sufficient thermalization
- Average over runs
- Test finite-size effects

## Community and Support
- Open source on GitHub
- Developed at UW-Madison (Prof. Jiamian Hu)
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/jhu238/AtomMag
2. Related publications from UW-Madison

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Limited
- Active development: Research code
- Specialized strength: GPU-parallel atomistic spin dynamics with 66x speedup
