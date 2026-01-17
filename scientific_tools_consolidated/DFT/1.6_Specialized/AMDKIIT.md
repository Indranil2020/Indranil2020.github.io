# AMDKIIT

## Official Resources
- Homepage: https://github.com/AMDKIIT/amdkiit
- Source Repository: https://github.com/AMDKIIT/amdkiit
- License: GNU General Public License v3.0

## Overview
AMDKIIT (`ab initio` Molecular Dynamics at KIIT/IIT) is a specialized Plane-Wave DFT software package developed to perform efficient molecular dynamics simulations. Created under the Indian National Supercomputing Mission (NSM), it is designed to run efficiently on high-performance computing clusters, including those with GPU acceleration. It bridges the gap between general-purpose DFT codes and specialized MD engines, focusing on the high-throughput generation of AIMD trajectories.

**Scientific domain**: Plane-Wave DFT, Ab Initio Molecular Dynamics
**Target user community**: HPC users, Researchers in materials chemistry and dynamical processes

## Theoretical Methods
- **Kohn-Sham DFT**: Standard formulations (LDA, PBE).
- **Plane-Wave Basis**: Periodic boundary conditions.
- **Pseudopotentials**: Norm-Conserving and Ultrasoft (UPF format support).
- **Car-Parrinello & Born-Oppenheimer MD**: Propagation schemes.
- **Thermostats**: Nose-Hoover, Berendsen for NVT ensembles.

## Capabilities (CRITICAL)
- **Electronic Minimization**: Self-consistent field using iterative diagonalization (Davidson/Conjugate Gradient).
- **Forces & Stress**: Analytic calculation of Hellmann-Feynman forces and stress tensors.
- **Molecular Dynamics**: Long-time scale NVE and NVT simulations.
- **GPU Acceleration**: Offloading of heavy FFT and BLAS operations to GPUs (CUDA).
- **Parallelism**: Hybrid MPI/OpenMP parallelization.

## Key Strengths

### GPU Optimization:
- Built from the ground up to leverage modern heterogeneous architectures (CPU+GPU).
- improved performance-per-watt for long MD runs.

### Local Development:
- Major indigenous code development project from India (IIT Kanpur).
- Open architecture allowing for academic contributions.

## Inputs & Outputs
- **Inputs**:
  - `input.in`: Main control file (cutoffs, convergence, MD steps).
  - `structure.xyz`: Initial coordinates.
  - Pseudopotentials (`.upf`).
- **Outputs**:
  - `trajectory.xyz`: MD Steps.
  - `energy.dat`: Thermodynamic logs.
  - `forces.dat`: Atomic forces.
  - `restart.bin`: Binary checkpoint capability.

## Interfaces & Ecosystem
- **File Formats**: Compatible with standard UPF pseudopotentials (Quantum ESPRESSO ecosystem).
- **Visualisation**: Trajectories readable by VMD, Ovito.

## Advanced Features
- **Berry Phase**: (Developmental) Polarization calculations.
- **Metadynamics**: (Planned) Enhanced sampling integration.

## Performance Characteristics
- **Speed**: Competitive with major codes for standard MD benchmarks on GPU nodes.
- **Scaling**: Good strong scaling on cluster partitions.

## Computational Cost
- **High Efficiency**: Design goal is to reduce wall-time for 10-100 ps simulations.

## Limitations & Known Constraints
- **Feature Set**: Less feature-rich than VASP/QE (e.g., no hybrid functionals, no GW yet).
- **Documentation**: Documentation is evolving; usage requires familiarity with standard PW-DFT inputs.
- **Maturity**: Newer code compared to established giants; expect rapid changes.

## Comparison with Other Codes
- **vs Quantum ESPRESSO**: Both use UPF/Plane-Waves; AMDKIIT is simpler and optimized specifically for MD on specific hardware.
- **vs CPMD**: CPMD is the ancestor of AIMD; AMDKIIT is a modern C++/CUDA implementation.
- **vs VASP**: VASP is the industry standard; AMDKIIT offers an open-source, GPU-ready alternative for basic MD tasks.
- **Unique strength**: GPU-native design philosophy for AIMD.

## Application Areas
- **Liquids & Solvation**: Structure of water, ions in solution.
- **Diffusion**: Ion migration in battery materials.
- **Surface Dynamics**: Adsorption and reconstruction processes.

## Best Practices
- **GPUs**: Running on CPU-only nodes misses the main optimization point.
- **Potentials**: Use standard GBRV or SSSP pseudopotentials (verify compatibility).
- **Time Step**: Use appropriate time steps (0.5 - 1.0 fs) for stability.

## Community and Support
- **Source**: Developed at IIT Kanpur.
- **GitHub**: Issues tracking available.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/AMDKIIT/amdkiit
2. NSM Project documentation.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GPLv3)
- Origin: Verified academic project.
