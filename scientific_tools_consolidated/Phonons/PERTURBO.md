# PERTURBO

## Official Resources
- Homepage: https://perturbo-code.github.io/
- Documentation: https://perturbo-code.github.io/mydoc_overview.html
- Source Repository: https://github.com/perturbo-code/perturbo
- License: GNU General Public License v3.0

## Overview
PERTURBO is an open-source software for first-principles calculations of charge transport and ultrafast carrier dynamics in materials with electron-phonon interactions. Developed at Caltech, PERTURBO computes electronic transport properties, carrier relaxation, and nonequilibrium dynamics using ab-initio electron-phonon matrix elements, handling spin-orbit coupling, polar materials, and providing comprehensive tools for studying electronic transport from first principles.

**Scientific domain**: Carrier transport, electron-phonon coupling, ultrafast dynamics  
**Target user community**: Transport properties, ultrafast spectroscopy, semiconductor physics

## Theoretical Methods
- Ab-initio electron-phonon coupling
- Boltzmann transport equation (BTE)
- Iterative solution of BTE
- Relaxation time approximation
- Carrier relaxation dynamics
- Ultrafast carrier dynamics
- Nonequilibrium distributions
- Temperature-dependent transport
- Spin-orbit coupling effects
- Polar corrections (Fröhlich)
- Wannier interpolation

## Capabilities (CRITICAL)
- Electronic transport coefficients from first principles
- Carrier mobility (electrons and holes)
- Electrical conductivity tensor
- Seebeck coefficient
- Electronic thermal conductivity
- Electron-phonon scattering rates
- Carrier relaxation times
- Ultrafast carrier dynamics
- Hot carrier cooling
- Nonequilibrium carrier distributions
- Temperature-dependent properties
- Spin-orbit coupling treatment
- Polar materials (Fröhlich interaction)
- Anisotropic transport tensors
- Integration with Quantum ESPRESSO
- HPC parallelization

**Sources**: Official PERTURBO documentation, Comp. Phys. Comm. 264, 107970 (2021)

## Key Strengths
- **Comprehensive transport**: Full iterative BTE solution beyond relaxation time approximation
- **Ultrafast dynamics**: Time-resolved carrier dynamics and hot carrier cooling
- **First-principles**: Ab-initio electron-phonon coupling, no empirical parameters
- **Modern implementation**: Efficient algorithms, HPC capable, Python interface

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO DFT output
  - Electron-phonon matrix elements
  - PERTURBO input files (pert.in)
  - Temperature and k-point lists
  
- **Output data types**:
  - Transport coefficients and mobility tensors
  - Scattering rates and relaxation times
  - Time-resolved carrier populations
  - Energy-resolved properties
  - Band-resolved contributions

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Native integration for DFT and phonon calculations
- **Wannier90**: Optional Wannier interpolation for enhanced efficiency
- **Python**: Python interface for post-processing and analysis
- **HDF5**: Efficient data storage for large calculations

## Workflow and Usage

### Typical Workflow:
```bash
# 1. DFT (Quantum ESPRESSO)
pw.x < scf.in > scf.out
pw.x < nscf.in > nscf.out

# 2. Phonon calculation
ph.x < ph.in > ph.out
ph.x < elph.in > elph.out

# 3. PERTURBO preprocessing
qe2pert.x -i qe2pert.in

# 4. PERTURBO calculation
perturbo.x -i pert.in > pert.out
```

### Mobility Calculation Input:
```
&perturbo
  calc_mode = 'trans'
  solver = 'iter'
  ftemper = 'silicon.temper'
  boltz_kdim = 40 40 40
  boltz_qdim = 20 20 20
  band_min = 1
  band_max = 8
/
```

## Advanced Features
- **Iterative BTE**: Full solution beyond RTA for accurate transport
- **Ultrafast dynamics**: Time-resolved populations and carrier thermalization
- **Polar materials**: Proper long-range Fröhlich interaction treatment
- **Spin-orbit coupling**: Spin-dependent scattering and transport

## Performance Characteristics
- **Computational cost**: DFT/phonon most expensive; PERTURBO efficient
- **Scalability**: HPC capable with MPI parallelization
- **k/q-grid**: Dense grids required for convergence
- **Typical runtime**: Hours to days depending on system and convergence

## Limitations & Known Constraints
- **Requires Quantum ESPRESSO**: DFT starting point necessary
- **Electron-phonon only**: Does not include electron-electron scattering
- **Convergence**: Multiple parameters require careful testing
- **Learning curve**: Moderate; requires understanding of transport theory
- **Platform**: Linux/Unix systems

## Comparison with Other Codes
- **vs EPW**: PERTURBO focuses on dynamics; EPW on superconductivity
- **vs BoltzTraP**: PERTURBO includes explicit electron-phonon; BoltzTraP uses constant τ
- **Unique strength**: Ultrafast carrier dynamics from first principles

## Application Areas
- **Semiconductors**: Carrier mobility, transport, device physics
- **Ultrafast spectroscopy**: Hot carrier dynamics, pump-probe theory
- **Optoelectronics**: Solar cells, LEDs, photodetectors
- **Thermoelectrics**: Transport coefficients, figure of merit

## Best Practices
- Quality DFT convergence and dense k/q-point grids
- Systematic convergence testing of all parameters
- Phonon calculation convergence critical
- Appropriate energy windows and band selection

## Community and Support
- Open-source (GPL v3)
- GitHub repository with active development
- Documentation website and user forum
- Workshop materials and tutorials

## Development
- Jin-Jian Zhou (lead developer, Caltech/IOP CAS)
- Marco Bernardi group (Caltech)
- Active development with regular updates

## Research Impact
PERTURBO enables first-principles calculations of carrier transport and ultrafast dynamics, advancing understanding of electronic transport and carrier relaxation in materials from ab-initio theory.

## Verification & Sources
**Primary sources**:
1. Homepage: https://perturbo-code.github.io/
2. Documentation: https://perturbo-code.github.io/mydoc_overview.html
3. GitHub: https://github.com/perturbo-code/perturbo
4. Publication: Comp. Phys. Comm. 264, 107970 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE and COMPREHENSIVE
- Documentation: DETAILED
- Source: OPEN (GitHub, GPL v3)
- Development: ACTIVE (Caltech)
- Applications: Ab-initio carrier transport, ultrafast dynamics, iterative BTE, electron-phonon scattering, hot carrier cooling, production quality
