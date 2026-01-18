# WOPTIC

## Official Resources
- **Homepage**: https://github.com/woptic/woptic
- **Documentation**: https://github.com/woptic/woptic/wiki
- **Source Repository**: https://github.com/woptic/woptic
- **License**: GNU General Public License (GPL)

## Overview
**WOPTIC** is a code designed to calculate frequency-dependent optical conductivity ($\sigma(\omega)$) and related transport properties using **Maximally Localized Wannier Functions (MLWFs)**. Its distinctive feature is the use of an **adaptive k-mesh refinement scheme** based on the tetrahedron method, which allows it to efficiently resolve fine spectral features (such as those arising from band crossings or van Hove singularities) that would require prohibitively dense uniform grids. It can also incorporate self-energy effects, making it suitable for correlated systems when coupled with DMFT.

**Scientific domain**: Optical Properties, Transport Theory, Correlated Electron Systems
**Target user community**: Researchers using WIEN2k and Wannier90 to study optical response in metals and semimetals

## Theoretical Methods
- **Kubo-Greenwood Formula**: Calculates the linear optical response.
- **Peierls Approximation**: For intraband contributions.
- **Wannier Interpolation**: Uses the tight-binding Hamiltonian in the Wannier basis to evaluate velocities and energies at arbitrary k-points.
- **Adaptive Mesh Refinement**: Iteratively subdivides tetrahedra in the Brillouin zone where the integrand (e.g., the spectral function) varies rapidly.
- **Many-Body Self-Energy**: Can include frequency-dependent self-energy $\Sigma(\omega)$ from DMFT calculations.

## Capabilities
- **Optical Conductivity**:
  - Real and Imaginary parts of $\sigma_{\alpha\beta}(\omega)$.
  - Decomposition into intraband and interband contributions.
- **DC Transport**:
  - DC Conductivity ($\sigma_{DC}$).
  - Seebeck Coefficient (Thermopower).
- **Plasma Frequency**: Calculation of Drude weights / plasma frequencies ($\omega_p$).
- **DMFT Integration**: Can read self-energy files to compute optical properties of correlated materials.

## Key Strengths
- **Convergence**: Achieves converged optical spectra with orders of magnitude fewer k-points than uniform grid methods.
- **Sharp Features**: Explicitly resolves Fermi surface sheets and band crossings.
- **Correlations**: One of the few public codes that combines Wannier interpolation with DMFT self-energies for optics.

## Inputs & Outputs
- **Inputs**:
  - Wannier Hamiltonian (`case_hr.dat`).
  - Position operator matrices (`case.mmn` or specific woptic format).
  - Dipole matrix elements.
  - Optional: Self-energy file (`sigma.dat`).
- **Outputs**:
  - Optical conductivity files (`opt_cond.dat`).
  - Plasma frequencies.
  - k-mesh refinement statistics.

## Interfaces & Ecosystem
- **Upstream**:
  - **WIEN2k**: Primary DFT engine.
  - **wien2wannier**: Interface to generate W90 inputs.
  - **Wannier90**: Generates the Wannier Hamiltonian.
- **Parallelism**: MPI parallelization for evaluating k-points.

## Performance Characteristics
- **Efficiency**: The adaptive scheme heavily reduces the computational cost for metals where the optical transitions are concentrated in specific regions of the BZ.
- **Scaling**: Scales linearly with the number of generated k-points (tetrahedra).

## Limitations & Known Constraints
- **Complexity**: Requires a working chain of WIEN2k -> wien2wannier -> Wannier90 -> WOPTIC.
- **Symmetry**: Does not fully exploit symmetry to reduce the irreducible BZ (calculates in full BZ usually), though adaptive mesh helps mitigate this.

## Comparison with Other Codes
- **vs. Wannier90 (postw90)**: Wannier90 also calculates optical conductivity (Berry phase module) but uses a uniform grid; WOPTIC's adaptive grid is superior for metals and sharp spectral features.
- **vs. TRIQS/DFTTools**: TRIQS also handles DMFT optics but WOPTIC's standalone adaptive mesh is a unique optimization for the integration.

## Application Areas
- **Iron Pnictides**: Optical response in correlated superconductors.
- **Semimetals**: Weyl and Dirac semimetals with sharp nodal points.
- **Transition Metals**: Transport properties of d-electron systems.

## Community and Support
- **Development**: Developed at TU Wien (Elias Assmann et al.).
- **Support**: Via GitHub issues.

## Verification & Sources
- **Official Website**: [https://github.com/woptic/woptic](https://github.com/woptic/woptic)
- **Primary Publication**: E. Assmann et al., Comput. Phys. Commun. 202, 1 (2016).
- **Verification status**: ✅ VERIFIED
  - Active repository.
  - Published methodology.
