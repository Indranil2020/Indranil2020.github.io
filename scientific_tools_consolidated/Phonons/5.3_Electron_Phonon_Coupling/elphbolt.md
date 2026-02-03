# elphbolt

## Official Resources
- **Repository**: https://github.com/nakib/elphbolt
- **License**: GPL-3.0

## Overview
**elphbolt** is a modern Fortran solver for the **coupled electron-phonon Boltzmann Transport Equations (BTE)**. While most transport codes assume one carrier system is in equilibrium (e.g., phonons are a thermal bath for electrons), elphbolt treats both electrons and phonons as non-equilibrium species that drag each other. This allows for the rigorous calculation of **phonon drag** contributions to the Seebeck coefficient and thermal transport in metals and semiconductors.

**Scientific domain**: Generalized Transport Theory, Thermoelectrics, Electron-Phonon Coupling
**Target user community**: Researchers studying phonon drag and hydrodynamic transport

## Theoretical Methods
- **Coupled BTE**: Simultaneously solves the linearized BTE for electrons ($f_{\mathbf{k}}$) and phonons ($n_{\mathbf{q}}$).
  - Electro-thermal driving forces ($\nabla T$, $\mathbf{E}$).
- **Scattering Matrices**:
  - Electron-Phonon (el-ph) scattering from ab initio (wannier interpolation).
  - Phonon-Phonon (ph-ph) scattering (3-phonon processes).
  - Electron-Impurity and Phonon-Isotope scattering.
- **Iterative Solver**: Self-consistent solution of the coupled distribution functions.

## Capabilities
- **Transport Coefficients**:
  - Electrical Conductivity ($\sigma$).
  - Seebeck Coefficient ($S$), including the explicit phonon drag component ($S_{drag}$).
  - Electronic Thermal Conductivity ($\kappa_e$).
  - Lattice Thermal Conductivity ($\kappa_{ph}$), modified by electron scattering.
- **Analysis**:
  - Spectral decomposition of drag.
  - Identification of hydrodynamic regimes.

## Key Strengths
- **Phonon Drag**: One of the very few public codes capable of calculating phonon drag from first principles, which is dominant in many semiconductors at low temperatures and in high-mobility materials.
- **Mutual Interactions**: Captures the effect of non-equilibrium electrons on phonon lifetimes (and vice versa).
- **Modern Fortran**: Written in object-oriented Fortran 2008 for modularity and performance.

## Inputs & Outputs
- **Inputs**:
  - Electron-phonon matrix elements (from EPW).
  - Phonon-phonon matrix elements (from Phono3py or similar).
  - Electronic/Phononic band structures (Wannier90).
- **Outputs**:
  - Transport coefficients vs Temperature.
  - Drag contributions.

## Interfaces & Ecosystem
- **EPW**: Designed to consume scattering data produced by the EPW code (Quantum ESPRESSO).
- **Wannier90**: Uses Wannier interpolation for dense grid sampling.


## Advanced Features
- **Coupled BTE**: Simultaneous electron and phonon non-equilibrium
- **Phonon drag**: Explicit drag contribution to Seebeck coefficient
- **Mutual scattering**: Electron-phonon and phonon-electron coupling
- **Hydrodynamic transport**: Momentum-conserving scattering regimes
- **Spectral decomposition**: Mode-resolved drag analysis
- **Modern Fortran**: Object-oriented Fortran 2008 implementation

## Computational Cost
- EPW input generation: Expensive (external)
- elphbolt BTE solution: Hours to days
- Memory: High (full scattering matrices)
- Scales with k-point and q-point grid density

## Best Practices
- Validate EPW electron-phonon matrix elements first
- Converge k-point and q-point grids systematically
- Check phonon-phonon scattering convergence
- Compare with experimental Seebeck coefficients
- Analyze temperature dependence of drag contribution

## Performance Characteristics
- **Computational Cost**: High. Dealing with full scattering matrices ($N_k \times N_q$) and coupling them requires significant memory and compute.
- **Parallelism**: MPI parallelization over k/q-points.

## Comparison with Other Codes
- **vs. EPW**: EPW calculates transport using the BTE for electrons (Ziman formula) but typically assumes equilibrium phonons (no drag, or simple models). elphbolt adds the full non-equilibrium phonon coupling.
- **vs. ShengBTE**: ShengBTE solves the phonon BTE (assuming equilibrium electrons). elphbolt couples the two.

## Application Areas
- **Thermoelectrics**: Determining the limit of ZT in materials where drag enhances thermopower.
- **Hydrodynamic Transport**: Studying electron fluids where momentum-conserving scattering dominates.
- **Pure Metals**: Low-temperature thermal conductivity anomalies.

## Community and Support
- **Development**: Caltech / Harvard (Nakib Protik).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/nakib/elphbolt](https://github.com/nakib/elphbolt)
- **Primary Publication**: N. H. Protik et al., Phys. Rev. B 102, 024302 (2020).
- **Verification status**: ✅ VERIFIED
  - State-of-the-art research code.
