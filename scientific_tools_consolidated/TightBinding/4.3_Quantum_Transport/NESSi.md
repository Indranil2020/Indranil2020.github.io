# NESSi

## Official Resources
- **Homepage**: https://nessi.tuxfamily.org/
- **Repository**: https://github.com/nessi-project (or SourceForge)
- **License**: GNU General Public License v3.0

## Overview
**NESSi** (Non-Equilibrium Systems Simulation) is an open-source C++ framework designed for the simulation of **non-equilibrium quantum dynamics** using the **Kadanoff-Baym** contour formalism. Unlike standard steady-state NEGF codes, NESSi explicitly works in the time domain, solving for the two-time Green's functions $G(t, t')$. This allows for the study of ultrafast phenomena, transient transport, and the relaxation dynamics of strongly interacting quantum systems driven by external fields.

**Scientific domain**: Ultrafast Dynamics, Many-Body Physics, Time-Dependent Transport
**Target user community**: Physicists studying non-equilibrium phases, optical control, and transient spectroscopy

## Theoretical Methods
- **Kadanoff-Baym Equations (KBE)**: Time propagation of the non-equilibrium Green's function on the L-shaped contour.
- **Many-Body Approximations**:
  - Hartree-Fock (HF).
  - Second Born Approximation (2B).
  - GW approximation.
  - T-Matrix approximations.
- **Interaction Types**:
  - Hubbard U (electron-electron).
  - Holstein model (electron-phonon).

## Capabilities
- **Simulations**:
  - Pump-probe spectroscopy experiments.
  - Quench dynamics.
  - Transient current in molecular junctions.
- **Observables**:
  - Time-resolved spectral functions $A(t, \omega)$.
  - Occupations and currents.
  - Energy conservation checks.
- **Systems**:
  - Lattice models (Hubbard, Holstein).
  - Single Anderson Impurity Model (SIAM).

## Key Strengths
- **Rigorous Dynamics**: Treats initial correlations and memory effects exactly within the chosen approximation (unlike GKBA which neglects some memory).
- **Modularity**: The software is structured as a library (`libnessi`), separating the Green's function data structures from the physical model implementation.
- **High-Order Solvers**: Advanced time-stepping algorithms for integro-differential equations.

## Inputs & Outputs
- **Inputs**:
  - Model parameters (hopping, U, coupling strength).
  - Time grid and contour definitions.
  - Initial state configuration.
- **Outputs**:
  - Two-time Green's functions (binary/HDF5).
  - Post-processed spectral data.

## Interfaces & Ecosystem
- **Libraries**: Depends on GSL, FFTW, MPI.
- **Usage**: Typically used to write specific solver executables for a given model Hamiltonian.

## Performance Characteristics
- **Computational Cost**: Scaling is $O(N_t^3)$ due to the two-time nature (memory effects), making it much more expensive than time-local methods (GKBA or TD-DFT).
- **Parallelism**: MPI parallelization over momentum/k-points or orbital blocks.

## Limitations & Known Constraints
- **Time Limits**: The cubic scaling limits simulations to relatively short physical times.
- **System Size**: Primarily suited for model Hamiltonians or small clusters; atomistic full-basis simulations are extremely costly in full KBE.

## Comparison with Other Codes
- **vs. TD-DFT**: TD-DFT is cheaper ($O(N_t)$) but often lacks memory effects and sophisticated correlation; NESSi captures lifetime effects and non-Markovian dynamics.
- **vs. CHEERS**: Another KBE solver; NESSi focuses on providing a C++ library infrastructure.

## Community and Support
- **Development**: Developed by Michael Schüler, Denis Golez, and collaborators (Graz/Fribourg).
- **Documentation**: Extensive Wiki and examples.

## Verification & Sources
- **Repository**: [https://nessi.tuxfamily.org](https://nessi.tuxfamily.org)
- **Primary Publication**: M. Schüler et al., Comp. Phys. Comm. 257, 107484 (2020).
- **Verification status**: ✅ VERIFIED
  - Published and documented research code.
