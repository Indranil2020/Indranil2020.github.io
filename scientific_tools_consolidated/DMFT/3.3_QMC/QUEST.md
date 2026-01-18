# QUEST (Quantum Electron Simulation Toolbox)

## Official Resources
- Homepage: http://physics.wm.edu/~shiwei/quest/ (Research Group)
- Documentation: Available within the package/group resources
- Source Repository: Contact Shiwei Zhang group (William & Mary / Flatiron CCQ)
- License: Academic/Research (typically on request)

## Overview
QUEST (Quantum Electron Simulation Toolbox) is a Determinant Quantum Monte Carlo (DQMC) code developed by Shiwei Zhang's group (College of William & Mary, now also Flatiron CCQ). It is designed for high-precision simulations of strongly correlated electron systems, particularly the 2D Hubbard model and its variants. QUEST implements efficient algorithms for ground-state and finite-temperature auxiliary-field QMC (AFQMC), emphasizing numerical stability and control over the sign problem through constrained path methods.

**Scientific domain**: Determinant QMC, AFQMC, Hubbard models, strongly correlated electrons
**Target user community**: Condensed matter theorists, research groups studying Hubbard physics

## Theoretical Methods
- Determinant Quantum Monte Carlo (DQMC)
- Auxiliary-Field Quantum Monte Carlo (AFQMC)
- Constrained Path Monte Carlo (CPMC)
- Ground-state projection
- Finite-temperature grand canonical ensemble
- Hubbard-Stratonovich transformation
- Numerical stabilization techniques

## Capabilities (CRITICAL)
- **Hubbard Model Solver**: Specialized for 2D Hubbard models (single and multi-orbital)
- **Sign Problem Control**: Implements constrained path approximation to mitigate sign problem
- **Numerical Stability**: Advanced matrix decomposition and stabilization
- **Observables**: 
  - Energy and double occupancy
  - Magnetic correlations (spin-spin)
  - Superconducting pairing correlations
  - Density matrices
- **Lattices**: Square, rectangular, and custom geometries
- **Performance**: Optimized for large-scale production runs

**Sources**: Shiwei Zhang group publications, "The Quantum Electron Simulation Toolbox (QUEST)"

## Inputs & Outputs
- **Input formats**:
  - Fortran/C-style input parameters (lattice size, interaction U, temperature, etc.)
  - Simulation control parameters (trotter step, number of walkers)
  
- **Output data types**:
  - Green's functions
  - Correlation functions
  - Statistical averages and error estimates
  - Checkpoint files

## Interfaces & Ecosystem
- **Research Code**: Often used as a standalone engine or integrated into group-specific workflows.
- **Related Tools**: Precursor/related to CPMC-Lab (Matlab package) and newer CCQ AFQMC libraries.
- **Development**: actively maintained within the research group context.

## Limitations & Known Constraints
- **Availability**: Not a public open-source repository like GitHub; typically obtained via collaboration or request.
- **Documentation**: Research-grade; assumes familiarity with DQMC/AFQMC.
- **Scope**: Primarily focused on Hubbard-type lattice models.

## Performance Characteristics
- **Scalability**: Designed for large-scale production runs on HPC clusters.
- **Stability**: Emphasizes numerical control over the sign problem, allowing access to lower temperatures than standard implementations.
- **Efficiency**: Optimized for 2D lattices.

## Comparison with Other Codes
| Feature | QUEST | ALF |
| :--- | :--- | :--- |
| **Method** | DQMC / CP-AFQMC | DQMC |
| **Sign Problem** | Constrained Path mitigation | Standard stabilization |
| **Language** | Fortran | Fortran |
| **Focus** | Ground state / Low T Hubbard | Finite T Lattice Models |
| **Availability** | Research Group / Collaborative | Open / Registration |

## Verification & Sources
**Primary sources**:
1. QUEST Project Page: http://physics.wm.edu/~shiwei/quest/
2. Shiwei Zhang Publications (e.g., Phys. Rev. Lett., Phys. Rev. B)
3. "Constraint Path Monte Carlo" methodology papers

**Confidence**: VERIFIED - Leading research code in the field

**Verification status**: ✅ VERIFIED
- **Category**: Academic Research Code
- **Status**: Active (used in diverse publications)
- **Developer**: Shiwei Zhang Group
- **Specialized strength**: Controlled accuracy for Hubbard models via Constrained Path MC.
