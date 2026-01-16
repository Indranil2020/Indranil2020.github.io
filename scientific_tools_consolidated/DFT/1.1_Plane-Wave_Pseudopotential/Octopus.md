# Octopus

## Official Resources
- Homepage: https://octopus-code.org/
- Documentation: https://octopus-code.org/wiki/Manual
- Source Repository: https://gitlab.com/octopus-code/octopus
- License: GPLv2+

## Overview
Octopus is a scientific software package for density functional theory (DFT) and time-dependent density functional theory (TDDFT). It relies on a **Real-Space Grid** discretization of the Kohn-Sham equations (Finite Difference), rather than a basis set like plane waves or Gaussians. This makes it exceptionally suited for time-dependent phenomena, optical properties, and response functions in finite systems and periodic solids.

**Scientific domain**: TDDFT, real-space DFT, optical properties, electron dynamics
**Target user community**: Theoretical physicists, irregular geometry researchers, TDDFT specialists

## Theoretical Methods
- Real-Space Grid (Finite Difference)
- Time-Dependent DFT (TDDFT)
- Kohn-Sham DFT
- Pseudopotentials (Norm-conserving)
- Real-time electron dynamics
- Quantum Optimal Control Theory

## Capabilities
- **Time-Dependent**:
  - Real-time propagation of electronic wavefunctions
  - Linear and non-linear optical response
  - High-harmonic generation
- **Ground State**:
  - DFT total energy and forces
  - Geometry optimization
- **Systems**:
  - Finite systems (molecules, clusters) with arbitrary boundaries
  - Periodic systems (1D, 2D, 3D)
- **Advanced**:
  - Non-collinear magnetism
  - Spin-orbit coupling
  - Superconducting DFT

## Key Strengths
- **Real-Space flexibility**: No periodic image interaction artifacts for charged/finite systems.
- **Time-Domain**: World-leading capabilities in real-time electron dynamics.
- **Parallelization**: Efficient domain decomposition on real-space grids.
- **Multi-scale**: Integated with various Maxwell solvers.

## Computational Cost
- **Scaling**: Good parallel scaling due to domain decomposition ($O(N)$ to $O(N^2)$ depending on operation).
- **Grid Density**: Cost increases rapidly with finer grid spacing (needed for hard potentials).
- **Vacuum**: Vacuum is cheaper than in PW codes (adaptive mesh or cutoffs), but standard grid fills volume.

## Best Practices
- **Grid Spacing**: Convergence is controlled by grid spacing (not energy cutoff); careful testing required.
- **Pseudopotentials**: Use Norm-Conserving potentials (HGH, SG15); Ultrasoft/PAW not standard.
- **Time-Step**: For TDDFT, time-step must be stable; use enforced symmetry for long propagations.

## Comparison with Other Codes
- **vs Plane-Wave (VASP/QE)**: Octopus is far superior for finite systems in strong fields or time-dependent dynamics. PW codes are better/faster for standard bulk ground states.
- **vs PARSEC**: Both are real-space; Octopus strongly emphasizes TDDFT/Dynamics, PARSEC emphasizes static electronic structure.

## Community and Support
- **Forum**: `octopus-users` mailing list.
- **Events**: Octopus schools and hackathons.
- **Development**: Active GitLab repository.

## Verification & Sources
**Primary sources**:
1. Official Website: https://octopus-code.org/
2. Tancogne-Dejean et al., J. Chem. Phys. 152, 124119 (2020).
3. Andrade et al., Phys. Chem. Chem. Phys. 17, 31371 (2015).

**Confidence**: VERIFIED - Major community code.
