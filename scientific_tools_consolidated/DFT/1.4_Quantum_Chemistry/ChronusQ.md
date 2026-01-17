# ChronusQ

## Official Resources
- Homepage: https://urania.chem.washington.edu/chronusq/
- Documentation: https://urania.chem.washington.edu/chronusq/wiki/
- Source Repository: https://github.com/liresearchgroup/chronusq_public
- License: GNU General Public License v3.0

## Overview
ChronusQ (Chronus Quantum) is an open-source ab initio electronic structure software designed for addressing complex problems requiring consistent treatment of time dependence, relativistic effects, many-body correlation, electron-nuclear coupling, and spin. Written in modern C++ with MPI/OpenMP parallelism.

**Scientific domain**: Relativistic quantum chemistry, excited states, time-dependent methods, magnetic properties  
**Target user community**: Researchers studying heavy elements, relativistic effects, and time-dependent phenomena

## Theoretical Methods
- Hartree-Fock (RHF, UHF, GHF, X2C-HF)
- Density Functional Theory (RKS, UKS, GKS, X2C-KS)
- Two-component relativistic methods (X2C)
- Four-component Dirac-Coulomb
- Real-time time-dependent DFT (RT-TDDFT)
- Configuration Interaction (CI)
- Coupled Cluster (CC) methods
- Gauge-including atomic orbitals (GIAOs)
- Finite magnetic field calculations
- Multi-reference methods

## Capabilities (CRITICAL)
- Ground and excited state calculations
- Two-component and four-component relativistic Hamiltonians
- Exact two-component (X2C) transformations
- Real-time propagation for excited states
- Magnetic properties (NMR, EPR, magnetizability)
- Finite magnetic field electronic structure
- Generalized Hartree-Fock (complex orbitals)
- Non-collinear spin-DFT
- Parallel execution (MPI + OpenMP)
- Modern C++ implementation with TiledArray tensor engine

## Key Strengths

### Relativistic Capabilities:
- Full four-component Dirac-Coulomb
- Exact two-component (X2C) transformation
- Spin-orbit coupling
- Picture-change corrections
- Heavy element chemistry

### Time-Dependent Methods:
- Real-time TDDFT
- Non-perturbative dynamics
- Strong-field phenomena
- Electronic stopping power
- Absorption spectra from propagation

### Magnetic Field Handling:
- Gauge-including atomic orbitals
- Finite magnetic field DFT
- Uniform magnetic fields
- NMR shielding tensors
- Magnetizability

### Modern Software Design:
- C++17 standard
- TiledArray tensor library
- MPI/OpenMP hybrid parallelism
- Modular architecture
- Extensible framework

## Inputs & Outputs
- **Input formats**:
  - ChronusQ input files
  - XYZ coordinate files
  - Basis set specifications
  
- **Output data types**:
  - Total energies and gradients
  - Molecular orbitals
  - Properties (dipoles, multipoles)
  - Density matrices
  - Time-dependent observables

## Interfaces & Ecosystem
- **Basis sets**: Standard Gaussian basis sets
- **Libraries**: TiledArray, BLAS/LAPACK, LibXC
- **Visualization**: Standard molecular formats
- **External codes**: Interface capabilities

## Advanced Features

### X2C Relativistic Methods:
- One-step X2C transformation
- Picture-change corrected properties
- Spin-orbit DFT
- Efficient for heavy elements

### Real-Time Dynamics:
- Predictor-corrector propagation
- Absorption spectra
- Electron dynamics
- Non-linear phenomena

### Generalized Methods:
- Complex orbitals (GHF/GKS)
- Non-collinear spin
- Broken symmetry solutions
- Kramers-restricted methods

## Performance Characteristics
- **Speed**: Efficient with TiledArray tensors
- **Accuracy**: High-level relativistic methods
- **System size**: Medium-sized molecules
- **Memory**: Distributed memory capable
- **Parallelization**: Excellent MPI/OpenMP scaling

## Computational Cost
- **HF/DFT**: Comparable to other codes
- **X2C**: Small overhead over non-relativistic
- **Four-component**: Higher cost (factor of 8-16)
- **RT-TDDFT**: Depends on propagation time
- **Typical**: Heavy element calculations efficient

## Limitations & Known Constraints
- **Active development**: Some features still maturing
- **Community size**: Smaller user base
- **Documentation**: Growing but not exhaustive
- **Analytic gradients**: Limited scope
- **Large systems**: Best for medium-sized molecules
- **Learning curve**: Relativistic methods require expertise

## Comparison with Other Codes
- **vs DIRAC**: Both relativistic; ChronusQ more RT-TDDFT focused
- **vs ReSpect**: ChronusQ broader scope, ReSpect NMR specialist
- **vs BAGEL**: Both modern C++; different method focus
- **vs ORCA**: ChronusQ more specialized in relativistic/time-dependent
- **Unique strength**: Unified treatment of relativity, time-dependence, and magnetism

## Application Areas

### Heavy Element Chemistry:
- Actinide and lanthanide complexes
- Spin-orbit effects
- Relativistic corrections
- Heavy metal catalysis

### Magnetic Properties:
- NMR chemical shifts
- EPR parameters
- Magnetizability
- Finite field calculations

### Ultrafast Dynamics:
- Attosecond phenomena
- Strong-field interactions
- Electronic stopping
- Photoionization

## Best Practices

### Relativistic Calculations:
- Use X2C for efficiency
- Four-component for benchmarks
- Appropriate basis sets for heavy elements
- Picture-change corrections for properties

### Time-Dependent Calculations:
- Appropriate time step
- Sufficient propagation time
- Perturbation strength checks
- Convergence monitoring

## Community and Support
- Open-source GPL v3
- Active GitHub development
- Li Research Group (U. Washington)
- Academic publications and citations
- Growing user community

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/liresearchgroup/chronusq_public
2. Li Research Group: https://urania.chem.washington.edu/chronusq/
3. Williams-Young et al., WIREs Comput. Mol. Sci. (2020) - ChronusQ paper
4. arXiv preprints on developments

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Documentation: Available
- Active development: Yes (Li Research Group)
- Academic citations: Growing
