# StochasticGW

## Official Resources
- Homepage: https://stochasticgw.github.io/
- Source Repository: https://github.com/stochasticGW/stochasticGW
- License: GNU General Public License v3.0

## Overview
StochasticGW is a specialized software package designed to perform GW calculations (Green's function G and screened Coulomb interaction W) for very large systems. It utilizes stochastic techniques to scale linearly with the system size, enabling the calculation of quasiparticle energies for systems with thousands to tens of thousands of electrons, which are computationally inaccessible to deterministic GW codes.

**Scientific domain**: Many-body perturbation theory, GW approximation, large-scale materials science
**Target user community**: Researchers studying quasiparticle energies in large nanostructures, defects, and disordered systems

## Theoretical Methods
- Stochastic resolution of identity
- Real-time propagation of stochastic orbitals
- Stochastic sampling of the Hilbert space
- GW approximation (G0W0)
- Separable resolution of the identity
- Time-dependent Hartree methods

## Capabilities (CRITICAL)
- Quasiparticle energies (Occupied and Unoccupied)
- Band gaps of large systems
- Linear scaling with system size O(N)
- Systems with >10,000 electrons
- Defect energy levels
- Disordered systems

**Sources**: Official website, LBL publications

## Key Strengths

### Linear Scaling:
- Stochastic sampling avoids N^3 or N^4 bottlenecks
- Scales linearly O(N) with electrons
- Enables 10,000+ electron calculations

### Large System Access:
- Nanocrystals
- Surface defects
- Interfaces
- Amorphous materials

### Algorithms:
- Real-time propagation reduces frequency integration cost
- Stochastic sampling replaces explicit summation

## Inputs & Outputs
- **Input formats**:
  - DFT wavefunctions/energies (from Quantum ESPRESSO, PARATEC)
  - Input parameters
  
- **Output data types**:
  - Quasiparticle energies
  - Self-energy corrections
  - Statistical error bars

## Interfaces & Ecosystem
- **DFT Codes**: Quantum ESPRESSO (v6.x supported), PARATEC
- **Language**: C++ / Fortran
- **Parallelization**: MPI / OpenMP

## Advanced Features

### Tensor Contraction:
- Stochastic orbitals replace dense tensor contractions
- Statistical averaging converges to exact result

### Variance Reduction:
- Techniques to minimize stochastic noise
- Efficient sampling of key subspaces

## Performance Characteristics
- **Speed**: Linear scaling, massive speedup for large N
- **Accuracy**: Converges to deterministic GW with sampling ^-1/2
- **System size**: 1000 - 10,000+ atoms
- **Parallelization**: Near-linear scaling with cores

## Computational Cost
- **High**: Requires many stochastic samples
- **Memory**: Lower than deterministic GW
- **Bottleneck**: Number of stochastic orbitals (N_zeta)
- **Trade-off**: CPU time vs statistical error bar

## Limitations & Known Constraints
- **Stochastic Noise**: Results have error bars
- **Convergence**: Slow sqrt(N) convergence of error
- **DFT Code**: Tied to specific DFT input formats
- **Method**: Primarily G0W0 (one-shot GW)

## Comparison with Other Codes
- **vs BerkeleyGW**: StochasticGW scales better for massive systems, BerkeleyGW exact for smaller
- **vs Yambo**: StochasticGW specialized for linear scaling regime
- **Unique strength**: Unmatched capability for 10,000+ electron GW calculations

## Application Areas
- **Nanowires**: Band gaps of realistic wires
- **Defects**: Color centers in large supercells
- **2D Materials**: Moiré patterns in twisted bilayers
- **Polymers**: Long disordered chains

## Best Practices
- **Sampling**: Check convergence of error bar with number of stochastic orbitals
- **Validation**: Compare small system result with deterministic code
- **DFT**: Ensure converged DFT ground state
- **Core count**: Use large HPC allocation for sampling

## Community and Support
- Open-source GPL v3
- Developed by Rabani/Neuhauser/Baer groups
- GitHub repository
- Manual and examples available

## Verification & Sources
**Primary sources**:
1. Website: https://stochasticgw.github.io/
2. GitHub: https://github.com/stochasticGW/stochasticGW
3. V. Vlcek et al., J. Chem. Theory Comput. 13, 3718 (2017)

**Confidence**: VERIFIED - Active academic code

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (GPL)
- Specialized strength: Linear-scaling GW for massive systems
