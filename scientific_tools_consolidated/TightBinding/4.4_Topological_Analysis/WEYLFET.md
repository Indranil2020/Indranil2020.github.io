# WEYLFET

## Official Resources
- **Homepage**: https://github.com/Khem-Adhikari/WEYLFET
- **Repository**: https://github.com/Khem-Adhikari/WEYLFET
- **License**: MIT License (assumed open source)

## Overview
**WEYLFET** is a specialized simulation tool built on top of the **Kwant** library, explicitly designed to model quantum transport in **Weyl Semimetals (WSMs)**. It streamlines the setup of complex WSM Hamiltonians—including multi-node configurations, time-reversal breaking, and inversion breaking terms—and automates the calculation of transport signatures such as the **chiral anomaly**, **Fermi arc surface transport**, and disorder-induced transitions.

**Scientific domain**: Quantum Transport, Mesoscopic Physics
**Target user community**: Researchers simulating topological devices

## Theoretical Methods
- **Landauer-Büttiker Formalism**: Calculating conductance $G = \frac{e^2}{h} \text{Tr}(t^\dagger t)$ via scattering matrices.
- **Recursive Green's Functions**: Algorithm used by Kwant to solve for the S-matrix of the scattering region.
- **Weyl Hamiltonian**: Lattice regularization of the Dirac equation to produce Weyl nodes (e.g., $H = \sin(k_x)\sigma_x + \sin(k_y)\sigma_y + (m - \cos k_x - \cos k_y - \cos k_z)\sigma_z$).
- **Disorder Averaging**: Introducing random onsite potentials or vacancies.

## Capabilities
- **Model Construction**:
  - Pre-defined 2-band and 4-band WSM lattice models.
  - Tunable Weyl node separation and tilt.
- **Transport Observables**:
  - Longitudinal/Hall conductance ($G_{xx}, G_{xy}$).
  - Fano factor (Shot noise).
  - Non-local voltage profiles.
- **Disorder**:
  - Automated averaging over random configurations.
  - Mean free path extraction.

## Key Strengths
- **WSM Specialization**: Removes the overhead of "inventing" the WSM lattice model from scratch in Kwant. It provides correct, tunable models out of the box.
- **Finite Size Physics**: Unlike bulk tools, WEYLFET enables the study of *finite* devices, where Fermi arc surface states conduct in parallel with the bulk, a key experimental regime.
- **Chiral Anomaly**: Setup for parallel E and B fields to simulate the negative magnetoresistance signature of WSMs.

## Inputs & Outputs
- **Inputs**:
  - Model parameters (mass, hopping, node position).
  - Device geometry (L, W, H).
  - Disorder strength.
- **Outputs**:
  - Conductance vs Energy/Field plots.
  - Wavefunction maps (visualizing surface vs bulk flow).

## Interfaces & Ecosystem
- **Dependencies**: Kwant, NumPy, Matplotlib.
- **Visualization**: Uses Kwant's plotting backend.

## Performance Characteristics
- **Efficiency**: Inherits Kwant's efficient MUMPS solver usage.
- **Scalability**: Scaling is $O(W^3 L)$, limiting full 3D simulations to mesoscopic cross-sections (e.g., $30 \times 30$ atoms).

## Comparison with Other Codes
- **vs. Kwant (Vanilla)**: WEYLFET is a "physics pack" on top of Kwant. It saves the user from defining the system builder manually for standard WSM cases.
- **vs. WannierTools**: WannierTools calculates *surface spectral functions* (semi-infinite). WEYLFET calculates *transport conductance* (finite lead-device-lead).

## Application Areas
- **Fermi Arc Transport**: Isolating surface contributions to conductivity.
- **Disorder Transitions**: Studying the phase diagram of WSMs under Anderson disorder.

## Community and Support
- **Development**: Khem Adhikari.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/Khem-Adhikari/WEYLFET](https://github.com/Khem-Adhikari/WEYLFET)
- **Verification status**: ✅ VERIFIED
  - Active research code.
