# PROFESS

## Official Resources
- Homepage: https://profess.dev/
- Documentation: https://profess.dev/
- Source Repository: https://github.com/profess-dev/profess-ad
- License: Open Source (GPL/MPL)

## Overview
PROFESS (PRinceton Orbital-Free Electronic Structure Software) is a high-performance Orbital-Free Density Functional Theory (OF-DFT) code. By optimizing the electron density directly without using wavefunctions (orbitals), it achieves linear scaling O(N) with system size and a very small prefactor, enabling quantum-mechanical simulations of systems with hundreds of thousands of atoms, particularly for main-group liquid metals and plasmas.

**Scientific domain**: Orbital-Free Density Functional Theory (OF-DFT), Large-scale Molecular Dynamics
**Target user community**: Materials scientists studying liquid metals, warm dense matter, and large scale defects

## Theoretical Methods
- Orbital-Free Density Functional Theory (OF-DFT)
- Kinetic Energy Density Functionals (KEDF) (Thomas-Fermi, von Weizsäcker, Wang-Teter, etc.)
- Local pseudopotentials
- Geometry optimization
- Ab initio molecular dynamics (AIMD)
- Periodic boundary conditions

## Capabilities (CRITICAL)
- Simulation of systems with >100,000 atoms
- Linear scaling O(N) computational cost
- Molecular dynamics (NVE, NVT)
- Geometry optimization
- Ion-electron potential evaluation
- Stress tensor calculations
- Auto-differentiation (PROFESS-AD)

## Key Strengths

### Extreme Scalability:
- Linear scaling with system size
- Handles systems inaccessible to Kohn-Sham DFT
- Comparison with classical force fields in speed (orders of magnitude slower but QM accuracy)

### Orbital-Free Methodology:
- No orbital diagonalization required
- Direct minimization of energy functional
- Efficient for simple metals (Al, Mg, Li, etc.)

### Modern Architecture (PROFESS-AD):
- PyTorch-based backend
- Automatic differentiation for gradients
- Easy implementation of new functionals
- GPU acceleration support

## Inputs & Outputs
- **Input**:
  - Structure files (ion positions)
  - Functional specifications
  - Pseudopotential files (BLPS)
- **Output**:
  - Total energies
  - Forces
  - Stress tensors
  - Optimized geometries
  - MD trajectories

## Interfaces & Ecosystem
- **Python**: PROFESS-AD provides a Python interface
- **PyTorch**: Built on PyTorch for flexibility
- **Post-processing**: Compatible with standard MD analysis tools

## Advanced Features
- **Auto-differentiation**: Facilitates derivative-based optimization
- **Kinetic Functionals**: Library of non-local and semi-local KEDFs
- **Warm Dense Matter**: High-temperature simulations

## Performance Characteristics
- **Speed**: Orders of magnitude faster than KS-DFT for large systems
- **System size**: 10,000 - 1,000,000+ atoms
- **Accuracy**: Dependent on KEDF quality; excellent for main-group metals

## Computational Cost
- **Scaling**: O(N)
- **Memory**: Low capabilities required primarily for density grid
- **Efficiency**: Very high for metallic systems

## Limitations & Known Constraints
- **KEDF Accuracy**: Accuracy limited by kinetic energy functional approximation
- **Transition Metals**: Difficult to treat accurately due to lack of shell structure in OF-DFT
- **Pseudopotentials**: Requires local pseudopotentials
- **Documentation**: Varying between versions (PROFESS 3 vs PROFESS-AD)

## Comparison with Other Codes
- **vs VASP/QE**: PROFESS is Orbital-Free (faster, less general)
- **vs DFTpy**: Similar OF-DFT focus, PROFESS specializes in high-performance C++/PyTorch
- **vs Classical MD**: PROFESS provides electronic structure effects lacking in force fields
- **Unique strength**: Massive scale AIMD for liquid metals

## Application Areas
- **Liquid Metals**: Structure and dynamics of liquid Al, Mg, Li
- **Warm Dense Matter**: High pressure/temperature physics
- **Defects**: Large-scale reconstruction in metals
- **Grain Boundaries**: Metallic interfaces

## Verification & Sources
**Primary sources**:
1. Official website: https://profess.dev/
2. L. Hung et al., J. Chem. Theory Comput. 9, 2196 (2013)
3. M. Chen et al., Comput. Phys. Commun. 190, 225 (2015)
4. GitHub: https://github.com/profess-dev/profess-ad

**Confidence**: VERIFIED
- Status: Active open source project (Carter Group)
- Existence: Confirmed
