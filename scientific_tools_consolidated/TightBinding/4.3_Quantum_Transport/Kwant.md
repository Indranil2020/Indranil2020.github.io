# Kwant

## Official Resources
- **Homepage**: https://kwant-project.org/
- **Documentation**: https://kwant-project.org/doc/
- **Repository**: https://gitlab.kwant-project.org/kwant/kwant
- **License**: BSD 2-Clause License

## Overview
**Kwant** is a powerful, open-source Python package for numerical quantum transport calculations. It allows for the construction of tight-binding models with arbitrary shapes and dimensionality and the calculation of their transport properties using the **scattering matrix** formalism (Landauer-Büttiker). Kwant is widely considered the community standard for mesoscopic transport due to its flexibility, ease of use, and "Builder" pattern which decouples the physics (Hamiltonian) from the geometry.

**Scientific domain**: Mesoscopic Physics, Topological Matter, Quantum Transport
**Target user community**: Theorists and experimentalists simulating quantum devices (QPCs, Hall bars, nanowires)

## Theoretical Methods
- **Tight-Binding**: Discrete lattice models with arbitrary hopping terms.
- **Scattering Matrix ($S$)**: Calculated via the wave-function matching method (or recursive Green's functions) for open systems connected to semi-infinite leads.
- **Landauer Formalism**: Conductance $G = \frac{2e^2}{h} \sum_n T_n$.
- **Green's Functions**: Integration for local quantities like density of states (LDOS).

## Capabilities
- **System Construction**:
  - Arbitrary lattices (1D, 2D, 3D) and complex shapes.
  - Symmetries (translational, rotational) handling.
- **Observables**:
  - Conductance and Shot Noise.
  - S-matrix elements (transmission/reflection amplitudes).
  - Wavefunctions $\psi(\mathbf{r})$ in the scattering region.
  - Local currents and spin densities.
- **Physics**:
  - Quantum Hall Effect (magnetic Peierls phases).
  - Superconductivity (Bogoliubov-de Gennes).
  - Topological insulators and Majorana modes.
  - Spin-Orbit Coupling.

## Key Strengths
- **Flexibility**: The "Builder" interface is extremely intuitive: `syst[site] = potential`.
- **Performance**: While the interface is Python, the heavy lifting is done by efficient C/Cython cores and sparse linear algebra (MUMPS, UMFPACK).
- **Ecosystem**: Highly extensible (e.g., `Tkwant` for time-dependent transport) and integrates perfectly with the SciPy stack.

## Inputs & Outputs
- **Inputs**: Python scripts defining the lattice, shape functions, and Hamiltonian values.
- **Outputs**:
  - S-matrices (NumPy arrays).
  - Scalar fields (LDOS) mapped to sites.
  - Band structures of leads.

## Interfaces & Ecosystem
- **Visualization**: Built-in plotting using Matplotlib (`kwant.plot`).
- **Tkwant**: Extension for time-dependent transport.
- **Qsymm**: Symmetry analysis of Hamiltonians.

## Performance Characteristics
- **Scaling**: Efficient for 2D/3D systems using sparse direct solvers ($O(N^{1.5...2.0})$ typically).
- **Parallelism**: Comparisons over parameters (energy, field) are embarrassingly parallel (e.g., `kwant.parallel`).

## Comparison with Other Codes
- **vs. OMEN/NEMO5**: Kwant is for model Hamiltonians (physics concepts); OMEN/NEMO5 are for atomistic material simulations (device engineering).
- **vs. Quantica.jl**: Quantica is a Julia-based spiritual successor/alternative to Kwant with higher performance for creating Hamiltonians but a smaller ecosystem.

## Application Areas
- **Majorana Fermions**: Simulating signatures of topological superconductivity in nanowires.
- **Quantum Hall**: Edge state transport and interference in interferometers.
- **Graphene**: Veselago lensing and p-n junction transport.

## Community and Support
- **Development**: CEA Grenoble (Xavier Waintal) and TU Delft (Anton Akhmerov).
- **Source**: GitLab.

## Verification & Sources
- **Website**: [https://kwant-project.org/](https://kwant-project.org/)
- **Primary Publication**: C. W. Groth et al., New J. Phys. 16, 063065 (2014).
- **Verification status**: ✅ VERIFIED
  - The Gold Standard for mesoscopic quantum transport.
