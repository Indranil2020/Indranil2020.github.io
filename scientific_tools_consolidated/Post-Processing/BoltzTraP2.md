# BoltzTraP2

## Official Resources
- Homepage: https://www.imc.tuwien.ac.at/forschungsbereiche_theorie_und_simulation/software_packages/boltztrap2/
- Documentation: https://gitlab.com/sousaw/BoltzTraP2
- Source Repository: https://gitlab.com/sousaw/BoltzTraP2
- License: GNU General Public License v3.0

## Overview
BoltzTraP2 is the modern Python-based successor to the widely used BoltzTraP code. It calculates electronic transport properties (Seebeck coefficient, conductivity, thermal conductivity, Hall coefficient) using the Boltzmann transport equation within the relaxation time approximation. BoltzTraP2 improves upon the original by using a more stable interpolation method, offering a Python API, and providing better integration with modern DFT workflows (via ASE, pymatgen).

**Scientific domain**: Electronic transport, thermoelectrics, Boltzmann transport  
**Target user community**: Thermoelectric researchers, materials scientists, Python users

## Theoretical Methods
- Semi-classical Boltzmann transport equation (BTE)
- Constant relaxation time approximation (CRTA)
- Smoothed Fourier interpolation (improved stability)
- Fermi integrals
- Transport distribution function
- Onsager coefficients

## Capabilities (CRITICAL)
- Calculation of transport coefficients (Seebeck, σ/τ, κe/τ)
- Fermi surface visualization (VTK export)
- Analysis of band structure derivatives (velocities, effective masses)
- Temperature and chemical potential dependence
- Python API for scripting and high-throughput workflows
- Command-line interface similar to original BoltzTraP
- Integration with ASE atoms and calculator outputs

**Sources**: BoltzTraP2 documentation, Comp. Phys. Comm. 231, 140 (2018)

## Inputs & Outputs
- **Input formats**: DFT output (VASP vasprun.xml, QE xml, Wien2k, Abinit, etc.), interpolation parameters
- **Output data types**: Transport coefficients (JSON/csv), Fermi surfaces (.vtk), plots

## Interfaces & Ecosystem
- **ASE**: Compatible with ASE calculators
- **Pymatgen**: Parsing via pymatgen
- **VASP, QE, Wien2k, Abinit, CASTEP**: Supported via parsers
- **Python**: Full library access

## Workflow and Usage
1. Run DFT calculation (dense k-mesh).
2. `bt2 interpolate` : Generate interpolation model from DFT output.
3. `bt2 integrate` : Calculate transport coefficients over T/mu grid.
4. `bt2 plot` : Visualize results.

## Performance Characteristics
- Computationally efficient interpolation
- Integration step is fast
- Memory usage can be high for very dense grids

## Application Areas
- Thermoelectric materials discovery
- Electronic conductivity in metals/semiconductors
- Hall effect
- Fermi surface topology

## Community and Support
- Open-source (GPL v3)
- Developed by Georg Madsen group (TU Wien)
- Active GitLab repository

## Verification & Sources
**Primary sources**:
1. Homepage: https://gitlab.com/sousaw/BoltzTraP2
2. Publication: G. K. H. Madsen, J. Carrete, M. J. Verstraete, Comp. Phys. Comm. 231, 140 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitLab)
- Development: ACTIVE
- Applications: Modern Python BTE solver, thermoelectrics
