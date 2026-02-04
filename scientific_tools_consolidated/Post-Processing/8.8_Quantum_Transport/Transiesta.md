# Transiesta

## Official Resources
- Homepage: https://departments.icmab.es/leem/siesta/
- Documentation: https://docs.siesta-project.org/projects/siesta/en/latest/tutorials/transiesta/index.html
- Source Repository: https://gitlab.com/siesta-project/siesta
- License: GNU General Public License v3.0

## Overview
Transiesta is the quantum transport module of the SIESTA density functional theory code. It uses the Non-Equilibrium Green's Function (NEGF) method combined with DFT to calculate electron transport properties of nanoscale systems under finite bias voltage. It enables the simulation of current-voltage characteristics, transmission spectra, and local currents in molecular junctions, nanowires, and interfaces.

**Scientific domain**: Quantum transport, NEGF-DFT, molecular electronics  
**Target user community**: Device physicists, nanoscientists, SIESTA users

## Theoretical Methods
- Non-Equilibrium Green's Function (NEGF)
- Density Functional Theory (DFT) with local orbitals (LCAO)
- Keldysh formalism for open systems
- Landauer-Büttiker formalism for conductance
- Pseudopotentials (Norm-conserving)
- Spin-polarized transport

## Capabilities (CRITICAL)
- **Transport Calculation**: I-V curves, zero-bias conductance, transmission function T(E)
- **Finite Bias**: Self-consistent calculation under applied voltage
- **Analysis**: Projected density of states (PDOS) in open systems, scattering states
- **Inelastic Transport**: Inelastic Electron Tunneling Spectroscopy (IETS) (via Inelastica)
- **Spin Transport**: Spin-polarized currents, spin torque
- **Scalability**: Capable of handling large systems (thousands of atoms) due to O(N) basis

**Sources**: SIESTA/Transiesta documentation, Comp. Phys. Comm. 147, 71 (2002)

## Inputs & Outputs
- **Input formats**: Flexible FDF format (SIESTA input), Electrode Hamiltonians (.TSHS)
- **Output data types**: Transmission files (.AVTRANS), PDOS, current density, eigenvalues

## Interfaces & Ecosystem
- **SIESTA**: Fully integrated part of the SIESTA package
- **TBTrans**: Post-processing tool for calculating transmission from Transiesta calculations
- **Inelastica**: Third-party tool for IETS and phonon effects using Transiesta
- **Sisl**: Python toolbox for manipulating Transiesta/SIESTA files

## Workflow and Usage
1. Calculate Electrode: Run SIESTA for the bulk electrode to get `electrode.TSHS`.
2. Setup Scattering Region: Define geometry (Left Electrode - Device - Right Electrode).
3. Run Transiesta: Perform self-consistent NEGF calculation (`transiesta input.fdf`).
4. Post-process: Run `tbtrans` to calculate transmission spectra and current.

## Performance Characteristics
- Highly efficient due to localized basis sets (linear scaling algorithms)
- Parallelized with MPI
- Efficient contour integration for Green's functions

## Application Areas
- Single-molecule transistors
- Graphene and 2D material nanoribbons
- STM tip-surface interactions
- Metallic contacts and interfaces
- Spin valves

## Community and Support
- Large SIESTA user community
- Active mailing list
- Developed by SIESTA developers (ICMAB, DTU, etc.)

## Verification & Sources
**Primary sources**:
1. Homepage: https://departments.icmab.es/leem/siesta/
2. Publication: M. Brandbyge et al., Phys. Rev. B 65, 165401 (2002)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitLab)
- Development: ACTIVE
- Applications: NEGF transport, finite bias, molecular electronics
