# WanTiBEXOS

## Official Resources
- **Repository**: https://github.com/ac-dias/wantibexos
- **Documentation**: https://wantibexos.readthedocs.io/
- **License**: GNU General Public License v3.0

## Overview
**WanTiBEXOS** (Wannier-based Tight-Binding for Excitonic and Optoelectronic Structures) is a Fortran90 code designed to efficiently compute the electronic, optical, and excitonic properties of materials ranging from bulk solids to nanostructures (0D, 1D, 2D). It bypasses the high computational cost of full *ab initio* GW-BSE calculations by solving the Bethe-Salpeter Equation (BSE) within a Wannier-based tight-binding framework, using parameters derived from **Wannier90**. This allows for the study of excitonic effects in large systems that would otherwise be inaccessible.

**Scientific domain**: Excitonics, Optoelectronics, 2D Materials, Nanostructures
**Target user community**: Researchers studying optical properties and excitons in complex or large unit-cell systems

## Theoretical Methods
- **Wannier Tight-Binding**: Uses the Maximally Localized Wannier Function (MLWF) Hamiltonian as the basis for the electronic structure.
- **Bethe-Salpeter Equation (BSE)**: Solves the effective two-particle Schrödinger equation for electron-hole pairs.
- **Tamm-Dancoff Approximation (TDA)**: Neglects the coupling between resonant and anti-resonant transitions (usually valid for optical absorption).
- **Screened Coulomb Potentials**: Implements dimensionality-dependent screening models (e.g., Keldysh potential for 2D, truncated potentials for 0D/1D) to account for dielectric environments without full GW screening.
- **Berry Phase**: Calculates topological invariants and Berry curvature integration.

## Capabilities
- **Electronic Properties**:
  - Band structures and Density of States (DOS).
  - Berry curvature distribution.
  - Topological invariants (Chern numbers).
- **Optical Properties**:
  - Dielectric function $\epsilon(\omega)$ (Real and Imaginary parts).
  - Absorption coefficient.
  - Refractive index.
  - Joint Density of States (JDOS).
- **Excitonic Properties**:
  - Exciton binding energies.
  - Exciton wavefunctions (real-space visualization).
  - Radiative lifetimes of excitons.
  - Exciton band structure.

## Key Strengths
- **Efficiency**: Orders of magnitude faster than *ab initio* BSE, making it suitable for high-throughput screening or large supercells.
- **Dimensionality Support**: Specific Coulomb interaction kernels for 3D (bulk), 2D (monolayers), 1D (nanotubes), and 0D (dots) systems.
- **Topological Analysis**: Integrated tools for analyzing topological phases alongside optical properties.
- **Parallelization**: OpenMP threading for efficient multi-core execution.

## Inputs & Outputs
- **Inputs**:
  - `wannier90_hr.dat`: Wannier Hamiltonian file.
  - `input.dat`: Main control file (system geometry, k-grid, potential type).
- **Outputs**:
  - `bands.dat`: Electronic bands.
  - `dielectric_function.dat`: Optical response.
  - `exciton_energies.dat`: Eigenvalues of the BSE Hamiltonian.
  - `wavefunctions`: Files for plotting exciton probability distributions.

## Interfaces & Ecosystem
- **Upstream**:
  - **Wannier90**: The primary source of the tight-binding Hamiltonian.
  - **DFT Codes**: Any code compatible with Wannier90 (QE, VASP, Siesta) can generate the starting point.
- **Visualization**: Outputs tailored for standard plotting tools (Gnuplot, Python).

## Performance Characteristics
- **Speed**: The tight-binding basis significantly reduces the matrix sizes compared to plane-wave BSE.
- **Scaling**: Efficient for large k-point grids required to converge excitonic spectra.
- **Memory**: Matrix storage scales with the number of Wannier orbitals squared/fourth power (for Two-Particle Hamiltonian), which is much smaller than plane-wave basis sizes.

## Limitations & Known Constraints
- **Screening Approximation**: Relies on model dielectric functions or static screening parameters rather than the full frequency-dependent $W(\omega)$ calculated in GW.
- **Basis Accuracy**: Results are only as good as the Wannierization; poor Wannier functions yield unphysical bands.
- **Correlation**: Does not include self-energy corrections to the band gap natively (typically requires a "scissor operator" or pre-corrected TB parameters).

## Comparison with Other Codes
- **vs. Yambo/BerkeleyGW**: WanTiBEXOS is not *ab initio* (uses TB models) and uses model screening, whereas Yambo/BerkeleyGW calculate full GW-BSE. WanTiBEXOS is much faster but less predictive for absolute band gaps without adjustments.
- **vs. PyGWBSE**: Both use TB-BSE approaches; WanTiBEXOS has extensive support for different dimensionalities (0D-3D).

## Application Areas
- **Transition Metal Dichalcogenides (TMDs)**: Excitons in monolayers (2D).
- **Perovskites**: Optical absorption in solar cell materials.
- **Topological Insulators**: Interplay of topology and excitonic effects.
- **Quantum Dots**: Optical properties of 0D confined systems.

## Community and Support
- **Development**: Developed by A.C. Dias (University of Brasilia).
- **Documentation**: Hosted on ReadTheDocs.
- **Issues**: GitHub issue tracker.

## Verification & Sources
- **Official Website**: [https://github.com/ac-dias/wantibexos](https://github.com/ac-dias/wantibexos)
- **Primary Publication**: A.C. Dias et al. (Check repository for specific citation, often associated with *Phys. Rev. B* papers on TMD excitons).
- **Verification status**: ✅ VERIFIED
  - Active maintained repository.
  - Clear documentation of physical models.
