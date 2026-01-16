# CP2K

## Official Resources
- Homepage: https://www.cp2k.org/
- Documentation: https://manual.cp2k.org/
- Source Repository: https://github.com/cp2k/cp2k
- License: GNU General Public License v2.0 or later

## Overview
CP2K is a versatile quantum chemistry and solid-state physics software package performing atomistic simulations of solid state, liquid, molecular, periodic, material, crystal, and biological systems. It excels at molecular dynamics simulations using mixed Gaussian and plane waves (GPW) method, and is particularly strong for large-scale condensed phase simulations including ab initio molecular dynamics.

**Scientific domain**: Molecular dynamics, DFT, condensed phase simulations, materials science  
**Target user community**: Researchers studying dynamics of molecules, liquids, solids, and biological systems

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Hybrid functionals (PBE0, HSE06, B3LYP)
- Gaussian and Plane Waves (GPW) method
- Quickstep DFT engine
- Born-Oppenheimer molecular dynamics (BOMD)
- Car-Parrinello molecular dynamics (CPMD)
- Time-Dependent DFT (TDDFT)
- Quantum Monte Carlo (QMC) interface
- Classical force fields (CHARMM, AMBER, GROMOS)
- QM/MM methods
- Semi-empirical methods (DFTB)

## Capabilities (CRITICAL)
- Ground-state electronic structure calculations
- Geometry optimization and transition state searches
- Ab initio molecular dynamics (Born-Oppenheimer, Ehrenfest)
- Path integral molecular dynamics (PIMD)
- Metadynamics and free energy calculations
- Linear-scaling DFT for large systems
- Hybrid DFT calculations with excellent performance
- Excited-state calculations via TDDFT
- Real-time TDDFT for non-linear spectroscopy
- NMR and EPR calculations
- Vibrational spectroscopy (IR, Raman)
- QM/MM simulations for biochemical systems
- Classical MD with machine learning potentials
- Surface and interface calculations
- Periodic and non-periodic boundary conditions
- Explicit solvation models
- Nudged elastic band (NEB) for reaction pathways
- GPU acceleration for selected methods

**Sources**: Official CP2K manual, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Main input file (.inp) with hierarchical keyword structure
  - Coordinate files (XYZ, PDB, CIF)
  - Restart files (.restart)
  - Basis sets and pseudopotentials (GTH, Goedecker-Teter-Hutter)
  
- **Output data types**:
  - Main output (.out) with energies, forces, stresses
  - Trajectory files (XYZ, DCD, PDB)
  - Restart files for continuation
  - Cube files for densities and orbitals
  - Atomic forces and stress tensors
  - Vibrational frequencies and modes

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - calculator interface
  - pymatgen - structure I/O
  - AiiDA - workflow automation
  - i-PI - path integral MD interface
  - PLUMED - enhanced sampling
  - libint - exact exchange integrals
  
- **QM/MM interfaces**:
  - GROMACS coupling
  - AMBER interface
  - CHARMM force fields
  
- **Post-processing**:
  - Built-in analysis tools
  - VMD visualization compatibility
  - TRAVIS trajectory analyzer
  
- **Machine learning**:
  - Neural network potentials integration
  - On-the-fly training capabilities

## Limitations & Known Constraints
- **Mixed basis complexity**: GPW method requires understanding of basis set convergence for both Gaussian and plane-wave components
- **Memory intensive**: Hybrid functionals and MP2 require substantial memory
- **Input complexity**: Hierarchical input format has steep learning curve
- **Efficient pseudopotentials**: Use GTH potentials optimized for GPW.
- **K-point sampling**: Limited support for k-points; best for large supercells
- **Scaling**: Parallel efficiency varies by method; hybrid DFT scales well
- **Documentation**: Extensive but can be difficult to navigate for beginners
- **Convergence**: Some systems require careful tuning of SCF parameters
- **Periodic/non-periodic mixing**: Some methods restricted to specific boundary conditions

## Best Practices
- **Basis Sets**: Use the "MOLOPT" basis sets (molecularly optimized) for best convergence.
- **Cutoff**: The `CUTOFF` parameter refers to the plane-wave grid; usually 300-600 Ry is needed. This is much higher than PW-only codes because it handles the density, not wavefunctions.
- **OT vs Diagonalization**: Use Orbital Transformation (OT) for large insulating systems (linear scaling); use diagonalization for metals.
- **Filesystem**: CP2K creates many temporary files; ensure fast I/O or use scratch directories.

## Verification & Sources
**Primary sources**:
1. Official website: https://www.cp2k.org/
2. Manual: https://manual.cp2k.org/
3. GitHub repository: https://github.com/cp2k/cp2k
4. J. Chem. Phys. 152, 194103 (2020) - CP2K: An electronic structure and molecular dynamics software package
5. Comput. Phys. Commun. 167, 103 (2005) - Quickstep: Fast and accurate DFT

**Secondary sources**:
1. CP2K tutorials and workshops
2. ASE calculator documentation
3. AiiDA-CP2K plugin
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Very active (mailing list, Google groups)
- Academic citations: >1,500 (main papers)
- Active development: Regular releases
