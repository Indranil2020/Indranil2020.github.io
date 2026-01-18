# dynamics-w90

## Official Resources
- **Repository**: https://github.com/michaelschueler/dynamics-w90
- **License**: GNU General Public License v3.0

## Overview
**dynamics-w90** is a sophisticated Fortran package designed to simulate **non-equilibrium electron dynamics** in solids using realistic tight-binding Hamiltonians derived from **Wannier90**. It focuses on light-matter interactions, enabling the study of ultrafast phenomena such as high-harmonic generation (HHG), time-resolved photoemission (tr-ARPES), and transient band structure engineering. The code distinguishes itself by implementing a **gauge-invariant formulation** for coupling electromagnetic fields to the Wannier basis, ensuring physical accuracy even in truncated basis sets.

**Scientific domain**: Ultrafast Spectroscopy, Non-linear Optics, Quantum Materials
**Target user community**: Theorists and experimentalists working on pump-probe spectroscopy and attosecond physics

## Theoretical Methods
- **Time-Dependent Schrödinger Equation (TDSE)**: Propagates the electronic wavefunction or density matrix in time under external fields.
- **Peierls Substitution + Corrections**: Implements a rigorous gauge-invariant coupling of the vector potential $\mathbf{A}(t)$ to the tight-binding Hamiltonian, including non-local terms.
- **Orbital Angular Momentum (OAM)**: Includes modern position operator corrections to calculating OAM and magnetic dichroism.
- **Berry Physics**: Evaluates topological quantities like Berry curvature and spin texture continuously in k-space.

## Capabilities
- **Spectroscopy Simulation**:
  - **High-Harmonic Generation (HHG)**: Calculates the time-dependent current $J(t)$ and its spectrum.
  - **Linear/Non-linear Optics**: Optical conductivity and higher-order susceptibilities.
  - **tr-ARPES**: Simulates time-resolved angle-resolved photoemission spectra (planned/beta).
- **Analysis**:
  - **Spin/Orbital Texture**: Maps spin and OAM on the Fermi surface.
  - **Population Dynamics**: Tracks excitation and relaxation of carriers.
- **System Types**: Bulk crystals, 2D materials (TMDs, Graphene), Topological Insulators.

## Key Strengths
- **Realism**: Uses DFT-derived parameters ($ab initio$ accuracy) rather than toy models.
- **Gauge Invariance**: Solves the long-standing problem of unphysical gauge-dependence in truncated basis simulations of optical response.
- **Efficiency**: Highly optimized Fortran 2008 core for time-propagation.
- **Modern Features**: Support for non-collinear spin and Spin-Orbit Coupling (SOC).

## Inputs & Outputs
- **Inputs**:
  - `_tb.dat` or `_hr.dat`: Wannier90 Hamiltonian.
  - `params.nml`: Namelist controlling the laser pulse (field strength, frequency, envelope) and time-stepping.
- **Outputs**:
  - `current.dat`: Time-dependent current (for HHG).
  - `populations.dat`: Time-dependent band populations.
  - `snapshots`: Wavefunction/Density matrix at specific times.

## Interfaces & Ecosystem
- **Upstream**: **Wannier90** (v2.x/v3.x compatible).
- **Helper Scripts**: Python scripts included for generating inputs and plotting results.

## Performance Characteristics
- **Computational Cost**: Scales linearly with the number of k-points and time steps; much cheaper than TD-DFT.
- **Parallelization**: MPI parallelization over k-points allows scaling to clusters.

## Limitations & Known Constraints
- **Correlations**: Currently primarily a mean-field / independent particle picture (with dephasing models). Full many-body scattering (Boltzmann/GKBA) is in development.
- **Basis Size**: Limited by the number of Wannier functions; extremely large bases (100+ bands) become expensive for time-propagation.

## Comparison with Other Codes
- **vs. TD-DFT (Octopus, Salmonella)**: dynamics-w90 is a "model" approach using fixed basis functions, allowing for much longer time scales and larger systems than full real-space TD-DFT.
- **vs. Python TB codes (Kwant)**: dynamics-w90 is specialized for *optical* driving and gauge-invariant field coupling, which is non-trivial in general TB codes.
- **vs. SBE codes**: Similar to Semiconductor Bloch Equation solvers but handles arbitrary multi-band structures from DFT.

## Application Areas
- **Valleytronics**: Selective valley excitation in TMDs.
- **Floquet Engineering**: Modifying topological properties with periodic driving.
- **Petahertz Electronics**: Analyzing current generation on sub-cycle time scales.

## Community and Support
- **Development**: Maintained by the group of Michael Schüler (PSI / University of Fribourg).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/michaelschueler/dynamics-w90](https://github.com/michaelschueler/dynamics-w90)
- **Primary Publication**: *Phys. Rev. B* 103, 125423 (2021) (Gauge invariance).
- **Verification status**: ✅ VERIFIED
  - Active research code.
  - Methodology published in high-impact journals.
