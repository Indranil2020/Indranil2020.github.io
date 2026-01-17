# MaZe

## Official Resources
- Homepage: https://gitlab.e-cam2020.eu/esl/MaZe
- Documentation: https://gitlab.e-cam2020.eu/esl/MaZe
- Source Repository: https://gitlab.e-cam2020.eu/esl/MaZe (E-CAM Gitlab)
- License: MIT License

## Overview
MaZe is a specialized code for performing Orbital-Free Density Functional Theory Molecular Dynamics (OF-DFT-MD) using the Mass-Zero (MaZe) constrained molecular dynamics approach. It focuses on the adiabatic propagation of the electronic density, treating it as a dynamic variable with zero mass, which allows for efficient and stable time evolution of large-scale systems within the orbital-free framework.

**Scientific domain**: Orbital-Free DFT, Molecular Dynamics, Algorithm Development
**Target user community**: Researchers in computational physics and extended Lagrangian dynamics

## Theoretical Methods
- Orbital-Free Density Functional Theory (OF-DFT)
- Mass-Zero Constrained Dynamics
- Extended Lagrangian formulation
- Adiabatic density propagation
- Kinetic Energy Density Functionals (Thomas-Fermi, etc.)
- Local pseudopotentials

## Capabilities (CRITICAL)
- OF-DFT Molecular Dynamics
- Constant energy (NVE) ensembles
- Efficient density propagation
- Large-scale metallic systems
- Integration with standard MD workflows

## Key Strengths

### Mass-Zero Algorithm:
- Avoids minimization of energy at every time step
- Propagates density constraints
- Maintains Born-Oppenheimer surface adherence
- Improved computational efficiency for dynamics

### Implementation:
- High-Performance Computing (HPC) ready
- C/C++ implementation
- MPI parallelization

## Inputs & Outputs
- **Input**:
  - Simulation parameters
  - Ionic configurations
  - Pseudopotentials
- **Output**:
  - Trajectories
  - Conserved quantities (Energy)
  - Density snapshots

## Interfaces & Ecosystem
- **E-CAM**: Part of the E-CAM software library for HPC
- **Libraries**: Uses standard numerical libraries (FFTW, BLAS)

## Advanced Features
- **Curvilinear Coordinates**: Generalization of MaZe to curvilinear constraints
- **Stability**: Enhanced stability over standard Born-Oppenheimer MD in some regimes

## Performance Characteristics
- **Speed**: Efficient propagation step
- **Scaling**: O(N) due to Orbital-Free nature
- **Accuracy**: Consistent with OF-DFT models

## Computational Cost
- **Complexity**: Linear with system size
- **Overhead**: Low overhead per timestep

## Limitations & Known Constraints
- **Method Scope**: Specialized for MaZe algorithm testing and usage
- **Functionals**: Limited to available OF-DFT functionals
- **Documentation**: Developer-focused

## Comparison with Other Codes
- **vs PROFESS**: PROFESS is a general purpose OF-DFT suite; MaZe focuses on the specific dynamics algorithm
- **vs CP2K**: CP2K uses orbital-based minimization (OT); MaZe uses orbital-free constraints
- **Unique strength**: Mass-Zero constrained dynamics implementation for OF-DFT

## Application Areas
- **Liquid Metal Dynamics**: Efficient sampling of phase space
- **Algorithm Research**: Testing extended Lagrangian methods
- **Large Scale MD**: Thousand-atom simulations

## Verification & Sources
**Primary sources**:
1. Repository: https://gitlab.e-cam2020.eu/esl/MaZe
2. References in E-CAM documentation
3. A. M. P. et al., J. Chem. Phys. (Related methodology papers)

**Confidence**: VERIFIED
- Status: Available Open Source
- Project: E-CAM Center of Excellence
