# elphbolt

## Official Resources
- **Repository**: https://github.com/nakib/elphbolt
- **License**: GNU General Public License v3.0

## Overview
**elphbolt** is a state-of-the-art solver for the **coupled electron-phonon Boltzmann Transport Equations (BTE)**. Unlike standard transport codes that treat electrons and phonons separately or use the relaxation time approximation, elphbolt solves the full system iteratively to capture complex non-equilibrium phenomena such as **phonon drag** (the dragging of electrons by a non-equilibrium phonon flux) and electron drag. It utilizes **Wannier90** for efficient interpolation of electronic bands and electron-phonon matrix elements, enabling first-principles calculations of thermoelectric and hydrodynamic transport in real materials.

**Scientific domain**: Thermoelectrics, Hydrodynamic Transport, Electron-Phonon Physics
**Target user community**: Researchers designing high-efficiency thermoelectric materials or studying non-equilibrium carrier dynamics

## Theoretical Methods
- **Coupled BTE**: Simultaneously solves the Boltzmann equations for electron distribution $f_{nk}$ and phonon distribution $n_{q\nu}$.
- **Wannier Interpolation**: Uses Maximally Localized Wannier Functions (MLWFs) to interpolate electron eigenvalues and electron-phonon matrix elements to dense k/q-grids.
- **Iterative Solver**: Self-consistently updates the distributions to converge the drag terms.
- **Eliashberg Theory**: Includes a module (`superconda`) for calculating phonon-mediated superconducting properties.

## Capabilities
- **Transport Coefficients**:
  - Electrical Conductivity ($\sigma$).
  - Seebeck Coefficient ($S$) with full phonon drag contribution.
  - Electronic Thermal Conductivity ($\kappa_e$).
  - Lattice Thermal Conductivity ($\kappa_{ph}$) with electron-phonon scattering.
- **Microscopic Analysis**:
  - Mode-resolved contributions to conductivity and drag.
  - Visualization of non-equilibrium distributions.
  - Lifetimes/Mean-free-paths for both carriers.
- **Superconductivity**: Calculation of $T_c$ and gap functions (via `superconda`).

## Key Strengths
- **Drag Physics**: One of the few public codes explicitly designed to capture phonon drag, which is dominant in many semiconductors at low/intermediate temperatures.
- **Ab Initio Accuracy**: No empirical parameters; inputs come directly from DFT/DFPT (e.g., Quantum ESPRESSO).
- **Efficiency**: OpenMP/MPI hybrid parallelization allows scaling to large grids necessary for transport convergence.

## Inputs & Outputs
- **Inputs**:
  - Wannier90 files (`_hr.dat`, `_chk`).
  - EPW files (`.epmatwp`) or similar electron-phonon element representations.
  - `elphbolt.in`: Control inputs.
- **Outputs**:
  - `onsager_coefficients.dat`: Full tensor of transport coefficients.
  - `linewidths`: Lifetimes of electrons/phonons.
  - `spectral_function`: Mode-resolved transport data.

## Interfaces & Ecosystem
- **Upstream**:
  - **Quantum ESPRESSO**: Generates the ground state and phonons (PHonon).
  - **EPW**: Often used to generate the initial coarse grid electron-phonon vertices.
  - **Wannier90**: Provides the tight-binding basis.
- **Parallelism**: Efficient handling of large memory requirements for coupled matrices.

## Performance Characteristics
- **Computational Cost**: More expensive than constant-time approximation codes (like BoltzTraP) but essential for accuracy in high-mobility or high-drag interactions.
- **Scaling**: Scales well with core count; memory bandwidth is often the bottleneck due to large interpolation tables.

## Limitations & Known Constraints
- **Input Complexity**: Requires a complex chain of preceding calculations (DFT -> PHonon -> Wannier90 -> EPW setup).
- **Validity**: BTE assumes well-defined quasiparticles; breaks down in the bad-metal or strongly correlated hopping regime (though elphbolt focuses on band transport).

## Comparison with Other Codes
- **vs. EPW**: EPW also calculates transport but elphbolt specializes in the *coupled* solution (drag), whereas standard EPW transport often uses frozen phonons or RTA.
- **vs. BoltzTraP**: BoltzTraP uses constant relaxation time (CRTA); elphbolt calculates fully energy- and momentum-dependent lifetimes from first principles.
- **vs. Perturbo**: Perturbo is similar (BTE solver) but elphbolt has a strong historical focus on the coupled drag effects.

## Application Areas
- **Thermoelectrics**: Optimization of $zT$ by engineering phonon drag contributions.
- **Hydrodynamic Electron Flow**: Regimes where electron-electron or electron-phonon scattering conserves momentum.
- **Low-Temperature Transport**: Analyzing the "phonon peak" in thermal conductivity.

## Community and Support
- **Development**: Developed by Nakib Protik (Harvard / MIT / now industry).
- **Source**: GitHub.
- **Citations**: Protik et al., Phys. Rev. B 102, 064302 (2020).

## Verification & Sources
- **Repository**: [https://github.com/nakib/elphbolt](https://github.com/nakib/elphbolt)
- **Primary Publication**: N. H. Protik, B. Kozinsky, et al.
- **Verification status**: ✅ VERIFIED
  - Active research code.
  - Verified against experimental data for Silicon and other standards.
