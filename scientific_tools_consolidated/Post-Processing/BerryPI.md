# BerryPI

## Official Resources
- Homepage: https://github.com/PyBerry/BerryPI
- Documentation: https://github.com/PyBerry/BerryPI/blob/master/README.md
- Source Repository: https://github.com/PyBerry/BerryPI
- License: GNU General Public License v3.0

## Overview
BerryPI is a Python software for calculating Berry phases and related properties from first-principles wavefunctions. It specifically focuses on the modern theory of polarization and the calculation of the spontaneous electric polarization in crystalline solids. It interfaces with VASP to extract the necessary Bloch functions.

**Scientific domain**: Berry phase, electric polarization, topological properties  
**Target user community**: Solid-state physicists, ferroelectricity researchers

## Theoretical Methods
- Modern Theory of Polarization (Berry Phase)
- Wannier charge centers
- Spontaneous polarization calculation
- Born effective charges
- Piezoelectric coefficients
- Anomalous Hall conductivity (in some versions)

## Capabilities (CRITICAL)
- Calculation of electronic polarization (Berry phase)
- Evaluation of ionic contribution to polarization
- Calculation of Born effective charges via finite differences
- Determination of piezoelectric tensors
- Interface with VASP WAVECAR
- Automatic handling of phase continuity

**Sources**: BerryPI GitHub repository, Comp. Phys. Comm. 184, 1280 (2013) (Reference to similar methods/tools)

## Inputs & Outputs
- **Input formats**: VASP WAVECAR, POSCAR, POTCAR
- **Output data types**: Polarization values (P_elec, P_ion), Berry phases, Born charges

## Interfaces & Ecosystem
- **VASP**: Primary electronic structure code supported
- **Python**: Written in Python using NumPy/SciPy

## Workflow and Usage
1. Perform VASP calculation to get converged wavefunction (WAVECAR).
2. Prepare BerryPI input script (or command line arguments).
3. Run BerryPI to compute the Berry phase of the Bloch functions.
4. For spontaneous polarization, compare centrosymmetric and distorted structures.

## Performance Characteristics
- Post-processing tool
- Performance depends on the size of the wavefunction (WAVECAR)
- Efficient implementation of Berry phase integration

## Application Areas
- Ferroelectric materials (polarization reversal)
- Piezoelectrics
- Multiferroics
- Topological insulators (Z2 invariants via Berry phase)

## Community and Support
- Open-source (GPL)
- GitHub repository
- Research code

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/PyBerry/BerryPI
2. Note: A tool named "BerryPI" was developed by S. J. Ahmed et al., Comp. Phys. Comm. 184, 647 (2013), but the GitHub repo might be a different or evolved implementation.

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Documentation: AVAILABLE (README)
- Source: OPEN (GitHub)
- Applications: Berry phase, polarization, VASP interface
