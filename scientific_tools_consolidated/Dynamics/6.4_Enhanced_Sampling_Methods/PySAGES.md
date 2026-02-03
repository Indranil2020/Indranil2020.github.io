# PySAGES

## Official Resources
- Homepage: https://github.com/SSAGESLabs/PySAGES
- Documentation: https://pysages.readthedocs.io/
- Source Repository: https://github.com/SSAGESLabs/PySAGES
- License: MIT

## Overview
PySAGES (Python Suite for Advanced General Ensemble Simulations) is a Python implementation of SSAGES that provides GPU-accelerated enhanced sampling methods. It offers a user-friendly Python interface and leverages JAX for automatic differentiation and GPU acceleration.

**Scientific domain**: GPU-accelerated enhanced sampling, free energy calculations  
**Target user community**: Researchers needing GPU-accelerated enhanced sampling with Python

## Theoretical Methods
- Metadynamics
- Adaptive biasing force (ABF)
- Umbrella sampling
- Harmonic bias
- Spectral ABF
- Neural network collective variables

## Capabilities (CRITICAL)
- GPU-accelerated sampling
- JAX backend
- Multiple MD engine support
- Neural network CVs
- Python interface
- Automatic differentiation

## Key Strengths

### GPU Acceleration:
- JAX backend
- Fast CV evaluation
- Efficient on GPUs
- Automatic differentiation

### Python Interface:
- Easy to use
- Flexible CVs
- Neural network support
- Modern API

## Inputs & Outputs
- **Input formats**:
  - Python configuration
  - MD engine inputs
  
- **Output data types**:
  - Free energy profiles
  - CV trajectories
  - Bias data

## Interfaces & Ecosystem
- **HOOMD-blue**: Integration
- **OpenMM**: Integration
- **LAMMPS**: Integration
- **JAX**: Backend

## Advanced Features
- **Neural network CVs**: Learned variables
- **Spectral ABF**: Improved ABF
- **GPU acceleration**: JAX-based
- **Autodiff**: Automatic gradients

## Performance Characteristics
- Excellent GPU performance
- Fast CV evaluation
- Good scaling
- Modern implementation

## Computational Cost
- GPU provides major speedup
- Efficient CV calculation
- Low overhead
- Overall: Excellent on GPU

## Best Practices
- Use GPU when available
- Validate CV choice
- Check convergence
- Use neural network CVs for complex systems

## Limitations & Known Constraints
- JAX dependency
- Newer than SSAGES
- Some methods still developing
- GPU recommended

## Application Areas
- Protein dynamics
- Chemical reactions
- Materials science
- Drug discovery

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/SSAGESLabs/PySAGES
2. P. Zubieta Rico et al., npj Comput. Mater. 10, 43 (2024)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
