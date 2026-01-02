# TurboRVB

## Official Resources
- Homepage: https://people.sissa.it/~sorella/web/index.html
- Documentation: https://turborvb.qe-forge.org/
- Source Repository: https://github.com/sissaschool/turborvb
- License: GNU General Public License v3.0

## Overview
TurboRVB is a quantum Monte Carlo package for ab initio electronic structure calculations with emphasis on correlated wavefunctions including Jastrow factors, geminals, and resonating valence bond (RVB) states. It provides advanced optimization techniques and efficient GPU implementation.

**Scientific domain**: Strongly correlated systems, high-accuracy QMC, GPU-accelerated calculations  
**Target user community**: Researchers studying strongly correlated materials, benchmark calculations

## Theoretical Methods
- Variational Monte Carlo (VMC)
- Lattice regularized diffusion Monte Carlo (LRDMC)
- Jastrow-Geminal wavefunctions
- Resonating valence bond (RVB) wavefunctions
- Antisymmetrized geminal power (AGP)
- Pfaffian wavefunctions
- Stochastic reconfiguration optimization
- Linear method optimization
- GPU-accelerated algorithms

## Capabilities (CRITICAL)
- Highly accurate ground-state energies
- Strongly correlated electron systems
- Advanced wavefunction ansatze (Jastrow-geminals, AGP)
- Efficient wavefunction optimization
- Structural optimization with QMC forces
- GPU acceleration for large systems
- Periodic boundary conditions
- Twisted boundary conditions
- Finite-size corrections
- Real-space electron densities
- Momentum distributions
- Superconducting order parameters
- Pseudopotential and all-electron calculations

**Sources**: Official TurboRVB documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - datasmin.input (main input)
  - fort.10 (wavefunction parameters)
  - Pseudopotentials
  - Initial trial wavefunctions
  
- **Output data types**:
  - out_min (optimization output)
  - fort.21 (energy output)
  - Optimized wavefunction parameters
  - Electron densities
  - Observable files

## Interfaces & Ecosystem
- **Wavefunction generation**:
  - Interface to DFT codes for initial guess
  - Quantum ESPRESSO interface
  
- **GPU support**:
  - CUDA implementation
  - Significant acceleration for large systems
  - Optimized for NVIDIA GPUs
  
- **Utilities**:
  - Conversion scripts
  - Analysis tools
  - Optimization utilities

## Limitations & Known Constraints
- **Specialized wavefunctions**: Requires understanding of advanced QMC ansatze
- **Computational cost**: QMC expensive; slower than DFT
- **GPU recommended**: CPU-only mode much slower
- **Statistical errors**: Stochastic method; results have error bars
- **Learning curve**: Steep; advanced QMC concepts required
- **Documentation**: Good but assumes QMC familiarity
- **Input complexity**: Many optimization parameters
- **System size**: GPU memory limits large systems
- **Platform**: Linux; GPU support requires CUDA

## Verification & Sources
**Primary sources**:
1. Official website: https://people.sissa.it/~sorella/web/
2. Documentation: https://turborvb.qe-forge.org/
3. GitHub repository: https://github.com/sissaschool/turborvb
4. K. Nakano et al., J. Chem. Phys. 152, 204121 (2020) - TurboRVB developments
5. S. Sorella et al., J. Chem. Phys. 143, 244112 (2015) - Geminal wavefunctions

**Secondary sources**:
1. TurboRVB tutorials and examples
2. Published strongly correlated system studies
3. GPU-accelerated QMC benchmarks
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (developers, SISSA)
- Academic citations: >200
- Unique strength: Advanced wavefunctions, GPU acceleration
