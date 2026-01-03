# EDRIXS

## Official Resources
- Homepage: https://github.com/NSLS-II/edrixs
- Documentation: https://nsls-ii.github.io/edrixs/
- Source Repository: https://github.com/NSLS-II/edrixs
- License: BSD 3-Clause License

## Overview
EDRIXS is an open-source toolkit for simulating Resonant Inelastic X-ray Scattering (RIXS) and X-ray Absorption Spectroscopy (XAS) based on exact diagonalization (ED) of model Hamiltonians. Developed at Brookhaven National Laboratory (NSLS-II), it allows users to model core-level spectroscopy of strongly correlated materials using crystal field theory and multiplet ligand field theory.

**Scientific domain**: X-ray spectroscopy (RIXS, XAS), exact diagonalization, strongly correlated systems  
**Target user community**: Spectroscopists, condensed matter physicists

## Capabilities (CRITICAL)
- **Exact Diagonalization**: Solves many-body Hamiltonians for small clusters/atoms.
- **Spectroscopy**: Calculates XAS, RIXS, and XPS spectra.
- **Hamiltonians**: Supports user-defined Hamiltonians involving Coulomb interactions, spin-orbit coupling, and crystal fields.
- **Geometry**: Handles experimental geometry (scattering angles, polarization).
- **Parallelization**: MPI-parallelized solver.

**Sources**: EDRIXS documentation, Comp. Phys. Comm. 241, 146 (2019)

## Inputs & Outputs
- **Input formats**: Python scripts defining parameters (Fk, Gk integrals, crystal field)
- **Output data types**: Spectra (intensity vs energy), ground state properties

## Interfaces & Ecosystem
- **Python**: Main interface.
- **Fortran/C++**: Core solvers.
- **Isotropic**: Database of atomic parameters.

## Workflow and Usage
1. Define atomic shell (e.g., Ni 3d).
2. Set parameters (Slater integrals, SOC, Crystal field).
3. Build Hamiltonian matrix.
4. Diagonalize to get eigenstates.
5. Calculate transition amplitudes (RIXS/XAS formula).
6. Plot spectrum.

## Performance Characteristics
- ED scales exponentially with Hilbert space size.
- Suitable for single-site or small cluster models.

## Application Areas
- Transition metal oxides
- High-Tc superconductors (cuprates)
- Interpreting synchrotron data

## Community and Support
- Developed by Y. L. Wang et al. (BNL)
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/NSLS-II/edrixs
2. Publication: Y. L. Wang et al., Comp. Phys. Comm. 241, 146 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: RIXS, XAS, Exact Diagonalization
