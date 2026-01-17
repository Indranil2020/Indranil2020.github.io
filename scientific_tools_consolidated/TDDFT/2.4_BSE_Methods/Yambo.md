# Yambo

## Official Resources
- Homepage: http://www.yambo-code.eu/
- Documentation: http://www.yambo-code.eu/wiki/
- Source Repository: https://github.com/yambo-code/yambo
- License: GNU General Public License v2.0

## Overview
Yambo is an open-source code for calculating excited state properties of materials from first principles using many-body perturbation theory. It implements the GW approximation for quasiparticle corrections and the Bethe-Salpeter equation for optical properties, featuring user-friendly interfaces and comprehensive capabilities for studying electronic excitations in molecules and solids.

**Scientific domain**: Many-body perturbation theory, GW, BSE, optical properties, excited states  
**Target user community**: Researchers studying electronic excitations, quasiparticle properties, and optical spectra

## Theoretical Methods
- GW approximation (G₀W₀, evGW, qsGW)
- Bethe-Salpeter equation (BSE)
- Time-Dependent Hartree-Fock (TDHF)
- Time-Dependent DFT (TDDFT)
- Dynamical Berry phase
- Real-time propagation
- Non-equilibrium Green's function (NEGF)
- Full-frequency integration (Real-axis, Godby-Needs, Plasmon-Pole)
- Hartree-Fock exchange
- Hybrid functionals
- Spin-orbit coupling

## Capabilities (CRITICAL)
- GW quasiparticle energies and band structures
- Optical absorption spectra including excitonic effects (BSE)
- Electron energy loss spectroscopy (EELS)
- Optical conductivity and dielectric function (RPA/BSE)
- Real-time propagation (improved in v5.3)
- Non-linear optics (SHG, DFG, multipole approximation - v5.2+)
- GPU acceleration (CUDA/OpenACC via devxlib - v5.3)
- Spin-orbit coupling and magnetic systems
- Interfaces with Quantum ESPRESSO, ABINIT
- Massively parallel (MPI/OpenMP/GPU)

**Sources**: Official Yambo documentation (v5.3), cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - DFT outputs from interfaced codes
  - Yambo input files (parameter-based)
  - Database files from DFT calculations
  
- **Output data types**:
  - Quasiparticle energies
  - Optical spectra
  - Dielectric functions
  - Self-energies
  - Screening functions
  - BSE eigenvalues and eigenvectors

## Interfaces & Ecosystem
- **DFT interfaces**:
  - Quantum ESPRESSO (primary)
  - ABINIT
  - PWscf
  - ETSF format support
  
- **Post-processing**:
  - yambopy - Python interface
  - Analysis scripts and utilities
  - Plotting tools
  
- **Workflow integration**:
  - AiiDA-Yambo plugin
  - Can be scripted for automated calculations

## Limitations & Known Constraints
- **Computational cost**: GW and BSE very expensive; limited to ~100-200 atoms
- **Memory intensive**: Self-energy and screening matrices large
- **DFT dependency**: Requires converged DFT ground state from external code
- **k-point convergence**: Often requires dense k-meshes
- **Frequency grid**: Convergence testing needed
- **Parallelization**: Complex; requires understanding of distribution
- **Learning curve**: Many-body methods require theoretical background
- **Documentation**: Good but assumes GW/BSE knowledge
- **Input complexity**: Many parameters to converge
- **Platform**: Primarily Linux/Unix; HPC recommended

## Verification & Sources
**Primary sources**:
1. Official website: http://www.yambo-code.org/
2. Documentation: http://www.yambo-code.org/wiki/
3. GitHub repository: https://github.com/yambo-code/yambo
4. A. Marini et al., Comput. Phys. Commun. 180, 1392 (2009) - Yambo code
5. D. Sangalli et al., J. Phys.: Condens. Matter 31, 325902 (2019) - Yambo developments

**Secondary sources**:
1. Yambo tutorials and schools
2. yambopy documentation
3. Published GW/BSE applications
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Very active (forum, schools, workshops)
- Academic citations: >800 (main papers)
- Active development: Regular releases, new features
