# Nanodcal

## Official Resources
- Homepage: https://www.nanoacademic.com/nanodcal
- Documentation: https://docs.nanoacademic.com/nanodcal/
- Source Repository: Proprietary (Nanoacademic)
- License: Proprietary / Academic License Available

## Overview
Nanodcal is a state-of-the-art quantum transport simulation software based on the Non-Equilibrium Green's Function (NEGF) density functional theory (DFT). It is designed to simulate electron transport through nanostructures and devices from first principles. Developed by Nanoacademic Technologies, it handles zero-bias and finite-bias conditions for molecular electronics, spintronics, and nanoscale devices.

**Scientific domain**: Quantum transport, NEGF-DFT, molecular electronics, spintronics  
**Target user community**: Device physicists, electrical engineers, materials scientists

## Theoretical Methods
- Non-Equilibrium Green's Function (NEGF) formalism
- Density Functional Theory (DFT) with LCAO basis
- Keldysh formalism for finite bias
- Spin-Orbit Coupling (SOC)
- Phonon scattering (inelastic transport)
- AC transport (time-dependent)

## Capabilities (CRITICAL)
- **Transport Properties**: I-V curves, transmission spectra, conductance, shot noise
- **Electronic Structure**: Band structure, DOS, complex band structure
- **Spin**: Collinear and non-collinear spin transport, spin torque
- **Device Simulation**: Two-probe systems (source-drain), multi-probe systems
- **Analysis**: Scattering states, transmission eigenstates, local currents
- **Thermal**: Thermoelectric coefficients (Seebeck), phonon transport (via NanoPhonon)

**Sources**: Nanodcal website, Phys. Rev. B 63, 245407 (2001) (Methodology)

## Inputs & Outputs
- **Input formats**: Python/Matlab-based scripting interface
- **Output data types**: HDF5 data, text files for I-V, transmission, etc.

## Interfaces & Ecosystem
- **Device Studio**: Graphical interface for building devices and analyzing results
- **Python**: Scripting API for automation
- **Parallelization**: MPI/OpenMP hybrid

## Workflow and Usage
1. Build device structure (Left Lead - Scattering Region - Right Lead).
2. Perform SCF calculation for leads (bulk).
3. Perform NEGF-SCF calculation for the central region (open system).
4. Calculate transmission and current.
5. Analyze results using Device Studio.

## Performance Characteristics
- Optimized for large-scale transport calculations
- Efficient handling of semi-infinite leads
- Parallelized for clusters

## Application Areas
- Molecular junctions
- Magnetic Tunnel Junctions (MTJ)
- 2D material transistors (FETs)
- Quantum point contacts
- Photocurrents

## Community and Support
- Commercial software with professional support
- Academic licenses available
- Developed by Nanoacademic Technologies (McGill University spinoff)

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.nanoacademic.com/nanodcal
2. Publication: J. Taylor, H. Guo, and J. Wang, Phys. Rev. B 63, 245407 (2001)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: PROPRIETARY
- Development: ACTIVE (Nanoacademic)
- Applications: NEGF-DFT transport, device simulation, I-V curves
