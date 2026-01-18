# Zen

## Official Resources
- Homepage: https://github.com/zen-dev/zen
- Documentation: https://github.com/zen-dev/zen (Project README and Wiki)
- Source Repository: https://github.com/zen-dev/zen
- License: Open Source

## Overview
Zen is a comprehensive computational toolkit developed for the *ab initio* simulation of strongly correlated materials. It is designed to seamlessly integrate Density Functional Theory (DFT) with Dynamical Mean-Field Theory (DMFT). The framework is built with a Julia-based core (`ZenCore`) for high-level orchestration and a Fortran-based engine for computationally intensive DMFT solving. It operates by manipulating parameters and data exchanged through configuration files, often orchestrating external DFT codes and internal solvers.

**Scientific domain**: Strongly correlated systems, DFT+DMFT simulations, Material science
**Target user community**: Researchers investigating correlated electrons, transition metal oxides, and exotic quantum phases

## Theoretical Methods
- Density Functional Theory (DFT) interfaces (VASP, Quantum ESPRESSO)
- Dynamical Mean-Field Theory (DMFT)
- Charge self-consistent DFT+DMFT
- Impurity Solvers (CT-HYB, NORG)
- Maximum Entropy Method (MaxEnt) for analytic continuation
- Julia interop for flexible workflow control

## Capabilities (CRITICAL)
- **Ab initio DFT+DMFT**: Fully integrated workflow for realistic materials.
- **Impurity Solvers**: Supports Continuous-Time Hybridization Expansion (CT-HYB) and Numerical Renormalization Group (NORG).
- **Charge Self-Consistency**: Updates electron density based on DMFT corrections of the charge density matrix.
- **Spectral Functions**: Computes spectral properties and Density of States (DOS) using analytic continuation tools (`ACFlow`).
- **Thermodynamics**: Calculation of free energy and thermodynamic properties.
- **Versatile Interface**: Connects with VASP and Quantum ESPRESSO for the DFT part of the cycle.

## Key Features

### Hybrid Architecture:
- **Julia Core (ZenCore)**: Provides a modern, high-level interface for workflow management and scripting.
- **Fortran Engine**: Ensures high performance for the heavy numerical lifting of the DMFT cycle (impurity solving).

### Integrated Solvers:
- Built-in support for advanced impurity solvers including CT-HYB and NORG.
- Designed to handle multi-orbital systems efficiently.

### Analytic Continuation:
- Includes tools (like `ACFlow`) for analytically continuing imaginary-axis data to real frequencies.

## Inputs & Outputs
- **Input formats**:
  - Julia scripts or configuration files defining model parameters (interaction U, J, inverse temperature $\beta$).
  - DFT output files (e.g., `WAVECAR`/`CHGCAR` from VASP, or specific Hamiltonian dumps).
  - Parameter blocks similar to DCore (e.g., `[model]`, `[impurity_solver]`, `[control]`) are often used in configuration files.
- **Output data types**:
  - Self-energies $\Sigma(i\omega_n)$
  - Green's functions $G(i\omega_n)$ and $G(\tau)$
  - Spectral functions $A(\omega)$
  - Thermodynamic observables (Occupancy, Energy)

## Interfaces & Ecosystem
- **DFT Integration**: Interfaces with VASP and Quantum ESPRESSO.
- **Julia Ecosystem**: Leverages Julia's scientific computing libraries (e.g., linear algebra, I/O).
- **Solver Interface**: Modular design allows for plugging in different impurity solvers.

## Workflow and Usage
The typical workflow involves:
1.  Running a DFT calculation (VASP/QE) to generate the initial non-interacting Hamiltonian and local basis.
2.  Setting up the DMFT cycle in Zen via a configuration script.
3.  `ZenCore` controls the iterative process:
    *   Mapping the lattice problem to an impurity model.
    *   Solving the impurity problem (Fortran engine).
    *   Updating the Self-energy.
    *   Solving the lattice Dyson equation.
    *   Updating the charge density (for self-consistency).
4.  Post-processing for spectral functions.

## Performance Characteristics
- **Efficiency**: Fortran backend ensures efficient handling of matrix operations and Monte Carlo steps.
- **Flexibility**: Julia frontend allows for easy customization and rapid prototyping of new workflows.
- **Parallelization**: Likely supports MPI/OpenMP for the impurity solver stage.


## Comparison with Other Frameworks
- **vs solid_dmft**: solid_dmft wraps TRIQS (Python/C++); Zen is a standalone Julia/Fortran toolkit.
- **vs TRIQS**: Zen aims for an integrated "all-in-one" experience; TRIQS is a modular toolbox requiring assembly.
- **Unique strength**: NORG solver for zero-temperature calculations and efficient Julia workflow.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/zen-dev/zen
2. Project Documentation/README
3. arXiv preprints associated with the development team (e.g., implementations of NORG, Zen framework papers)

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Active development: Recent commits observed
- Integration: Interfaces with standard DFT codes
