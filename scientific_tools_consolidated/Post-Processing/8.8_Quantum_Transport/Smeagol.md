# Smeagol

## Official Resources
- Homepage: https://www.tcd.ie/Physics/Smeagol/
- Documentation: https://www.tcd.ie/Physics/Smeagol/manuals.html
- Source Repository: Distributed via website (Academic license)
- License: Academic License (Free)

## Overview
Smeagol is a software package for calculating spin-dependent electron transport in nanoscale devices. It interfaces with the SIESTA DFT code to perform Non-Equilibrium Green's Function (NEGF) calculations. Smeagol is specifically designed for spintronics, handling magnetic materials, spin torque, and non-collinear magnetism in transport junctions.

**Scientific domain**: Quantum transport, spintronics, NEGF-DFT  
**Target user community**: Spintronics researchers, magnetics community

## Theoretical Methods
- Non-Equilibrium Green's Function (NEGF)
- Density Functional Theory (via SIESTA)
- Spin-dependent transport
- Non-collinear magnetism
- Spin-Transfer Torque (STT)
- Inelastic scattering (weak electron-phonon coupling)

## Capabilities (CRITICAL)
- **Transport**: I-V characteristics, transmission coefficients, magnetoresistance
- **Spin**: Spin-polarized currents, giant magnetoresistance (GMR), tunneling magnetoresistance (TMR)
- **Torque**: Calculation of spin-transfer torque vectors
- **Magnetic**: Self-consistent calculation of magnetic moments under bias
- **Interfaces**: Metal-molecule-metal junctions, magnetic tunnel junctions

**Sources**: Smeagol website, Phys. Rev. B 73, 085414 (2006)

## Inputs & Outputs
- **Input formats**: SIESTA fdf files, electrode Hamiltonians
- **Output data types**: Current-voltage curves, transmission spectra, density matrices

## Interfaces & Ecosystem
- **SIESTA**: Requires SIESTA as the DFT engine
- **Plotting**: Standard text output for plotting

## Workflow and Usage
1. Calculate bulk electrodes using SIESTA.
2. Setup transport geometry (leads + scattering region).
3. Run Smeagol (modified SIESTA executable) to solve NEGF self-consistently.
4. Extract transmission and currents.

## Performance Characteristics
- Computationally intensive (NEGF inversion)
- Parallelized (MPI)
- Memory scaling depends on system cross-section

## Application Areas
- Molecular spintronics
- Magnetic Tunnel Junctions (MTJ)
- Spin filters
- Graphene spintronics

## Community and Support
- Developed by Sanvito Group (Trinity College Dublin)
- Academic user base
- Regular workshops (Computational Spintronics)

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.tcd.ie/Physics/Smeagol/
2. Publication: A. R. Rocha et al., Phys. Rev. B 73, 085414 (2006)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: ACADEMIC (Registration)
- Development: ACTIVE (Sanvito Group)
- Applications: Spintronics, NEGF, magnetic transport
