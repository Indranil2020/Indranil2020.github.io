# Phoebe

## Official Resources
- Homepage: https://mir-group.github.io/phoebe/
- Documentation: https://phoebe.readthedocs.io/
- Source Repository: https://github.com/mir-group/phoebe
- License: Apache License 2.0

## Overview
Phoebe is a modern, high-performance code for calculating phonon and electron thermal transport properties from first principles. Developed at MIT, Phoebe solves the Boltzmann transport equation for both phonons and electrons with electron-phonon coupling, focusing on computational efficiency and advanced transport phenomena. The code features GPU acceleration, advanced algorithms, and handles coupled electron-phonon transport in a unified framework.

**Scientific domain**: Thermal transport, thermoelectrics, coupled electron-phonon dynamics  
**Target user community**: Thermal transport researchers, thermoelectric materials, computational materials science

## Theoretical Methods
- Phonon Boltzmann transport equation (BTE)
- Electron Boltzmann transport equation
- Coupled electron-phonon transport
- Iterative and variational BTE solutions
- Relaxation time approximation
- Phonon-phonon scattering (3-phonon processes)
- Electron-phonon scattering
- Phonon-boundary scattering
- Phonon-isotope scattering
- Electron-impurity scattering
- Wannier interpolation for electrons

## Capabilities (CRITICAL)
- Lattice thermal conductivity from first principles
- Electronic thermal conductivity
- Coupled electron-phonon thermal transport
- Electrical conductivity
- Seebeck coefficient
- Phonon and electron lifetimes
- Spectral thermal conductivity
- Cumulative thermal conductivity
- Mode-resolved transport properties
- Temperature-dependent transport
- Nanostructure and boundary scattering
- GPU acceleration for large systems
- HPC parallelization (MPI + OpenMP)
- Wannier interpolation for efficient calculations

**Sources**: Official Phoebe documentation, Nature Communications 12, 2222 (2021)

## Key Strengths
- **GPU acceleration**: Significant speedup for large-scale calculations
- **Coupled transport**: Unified electron-phonon treatment
- **Modern architecture**: C++ with Python interface, HPC-optimized
- **Advanced algorithms**: Variational and iterative BTE solvers

## Inputs & Outputs
- **Input formats**:
  - Phonopy force constants (for phonons)
  - Wannier90 data (for electrons)
  - Electron-phonon matrix elements
  - Crystal structure files
  - Phoebe configuration files
  
- **Output data types**:
  - Thermal conductivity tensors
  - Transport coefficients
  - Scattering rates and lifetimes
  - Spectral and cumulative properties
  - Mode-resolved contributions

## Interfaces & Ecosystem
- **Phonopy**: Import harmonic phonon properties
- **Quantum ESPRESSO**: Via phonopy interface for phonons
- **Wannier90**: For electronic structure and electron-phonon coupling
- **Python**: Python interface for workflow automation
- **HDF5**: Efficient data storage and exchange

## Workflow and Usage

### Phonon Transport Workflow:
```bash
# 1. Prepare phonon force constants (from phonopy/phono3py)
# 2. Create Phoebe input
# 3. Run Phoebe
mpirun -np 16 phoebe -in input.phoebe

# GPU acceleration
phoebe -in input.phoebe --useGPU
```

### Coupled Electron-Phonon Transport:
```bash
# Requires both phonon and Wannier electron data
phoebe -in coupled_transport.phoebe
```

## Advanced Features
- **Variational BTE**: Beyond relaxation time approximation
- **GPU kernels**: Optimized CUDA kernels for scattering calculations
- **Adaptive grids**: Smart q-point and k-point sampling
- **Nanostructures**: Boundary and grain scattering models
- **Hydrodynamic phonons**: Advanced transport regimes

## Performance Characteristics
- **GPU speedup**: 10-100x faster than CPU-only for large systems
- **Parallelization**: Excellent scaling with MPI+OpenMP
- **Memory**: Optimized for large k/q grids
- **Typical runtime**: Hours with GPU; days CPU-only for production

## Computational Cost
- Force constant calculations (DFT) most expensive
- Phoebe very efficient with GPU acceleration
- Iterative BTE more expensive than RTA
- Dense grids feasible with GPU

## Limitations & Known Constraints
- **GPU recommended**: CPU-only slower for large systems
- **Requires force constants**: From external phonon codes
- **Learning curve**: Moderate; requires transport theory knowledge
- **Documentation**: Growing; some advanced features need expertise
- **Platform**: Linux; GPU support requires CUDA

## Comparison with Other Codes
- **vs ShengBTE/phono3py**: Phoebe has GPU acceleration and coupled transport
- **vs ALAMODE**: Phoebe focuses on transport solver efficiency
- **Unique strength**: GPU-accelerated coupled electron-phonon transport

## Application Areas
- **Thermoelectrics**: Figure of merit (ZT), optimizing transport properties
- **Thermal management**: Heat dissipation in electronics
- **Nanostructures**: Grain boundaries, interfaces, thin films
- **Novel materials**: Materials with complex transport physics
- **High-throughput**: Rapid screening with GPU acceleration

## Best Practices
- Use GPU acceleration for production calculations
- Converge k/q-point grids systematically
- Test RTA vs iterative BTE convergence
- Validate against experimental data when available
- Appropriate boundary scattering parameters for nanostructures

## Community and Support
- Open-source (Apache 2.0)
- GitHub repository
- Documentation website
- MIT development team
- Growing user community
- Research collaborations

## Educational Resources
- Comprehensive documentation
- Tutorial examples
- Publication describing methodology
- Example input files
- Python API examples

## Development
- MIT Materials Intelligence Research group
- Active development
- GPU optimization ongoing
- Feature additions for advanced transport
- Community contributions

## Research Impact
Phoebe enables efficient first-principles thermal transport calculations with GPU acceleration, particularly valuable for coupled electron-phonon transport in thermoelectric materials and high-throughput materials screening.

## Verification & Sources
**Primary sources**:
1. Homepage: https://mir-group.github.io/phoebe/
2. Documentation: https://phoebe.readthedocs.io/
3. GitHub: https://github.com/mir-group/phoebe
4. Publication: Nature Communications 12, 2222 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub, Apache 2.0)
- Development: ACTIVE (MIT)
- Applications: GPU-accelerated thermal transport, coupled electron-phonon BTE, variational transport solvers, thermoelectrics, high-performance computing, production quality
