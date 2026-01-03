# BoltzTraP

## Official Resources
- Homepage: https://www.imc.tuwien.ac.at/forschungsbereiche_theorie_und_simulation/software_packages/boltztrap/
- Documentation: Included in source distribution
- Source Repository: https://www.imc.tuwien.ac.at/forschungsbereiche_theorie_und_simulation/software_packages/boltztrap/
- License: GNU General Public License v3.0

## Overview
BoltzTraP (Boltzmann Transport Properties) is a code for calculating electronic transport properties (Seebeck coefficient, electrical conductivity, electronic thermal conductivity) as a function of temperature and chemical potential. It solves the semi-classical Boltzmann transport equation within the constant relaxation time approximation (RTA), using a smoothed Fourier interpolation of the bands.

**Scientific domain**: Electronic transport, thermoelectrics, Boltzmann transport equation  
**Target user community**: Thermoelectric materials researchers, solid-state physicists

## Theoretical Methods
- Semi-classical Boltzmann transport equation
- Constant relaxation time approximation (RTA)
- Rigid band approximation
- Smoothed Fourier interpolation of band energies
- Calculation of transport distribution function tensor

## Capabilities (CRITICAL)
- Calculation of Seebeck coefficient (thermopower)
- Electrical conductivity (relative to relaxation time τ)
- Electronic thermal conductivity
- Power factor
- Hall coefficient
- Temperature and doping dependence of transport properties
- Interfaces with VASP, WIEN2k, Quantum ESPRESSO, CASTEP, etc.

**Sources**: BoltzTraP website, Comp. Phys. Comm. 175, 713 (2006)

## Inputs & Outputs
- **Input formats**: `intra` (structure/energy window), `energy` (eigenvalues), `struct` (WIEN2k structure format)
- **Output data types**: `.trace` (transport vs μ/T), `.condtens` (tensors), `.sig` (conductivity spectral function)

## Interfaces & Ecosystem
- **WIEN2k**: Native interface (`x_trans BoltzTraP`)
- **VASP**: Can be used via `vasp2boltz.py` scripts
- **Quantum ESPRESSO**: Interface available
- **Boltztrap2**: Modern Python successor (separate tool)

## Workflow and Usage
1. Perform DFT calculation with dense k-mesh.
2. Prepare input files (`case.intrans`, `case.energy`).
3. Run `BoltzTraP`.
4. Analyze output files (`case.trace`).

## Performance Characteristics
- Very fast (minutes) compared to DFT
- Fourier interpolation step scales with number of k-points and bands

## Application Areas
- Thermoelectric materials screening
- Doping optimization
- Transport in metals and semiconductors
- Hall effect simulations

## Community and Support
- Developed by Georg Madsen (TU Wien) and David Singh
- Widely cited legacy code
- Superseded by BoltzTraP2 for new projects

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.imc.tuwien.ac.at/forschungsbereiche_theorie_und_simulation/software_packages/boltztrap/
2. Publication: G. K. H. Madsen and D. J. Singh, Comp. Phys. Comm. 175, 713 (2006)

**Confidence**: VERIFIED

**Verification status**: VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GPL)
- Development: STABLE (Legacy)
- Applications: Thermoelectrics, transport, RTA
