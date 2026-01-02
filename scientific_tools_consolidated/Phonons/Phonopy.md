# Phonopy

## Official Resources
- Homepage: https://phonopy.github.io/phonopy/
- Documentation: https://phonopy.github.io/phonopy/
- Source Repository: https://github.com/phonopy/phonopy
- License: BSD 3-Clause License

## Overview
Phonopy is the standard code for phonon calculations using the finite displacement method with forces from DFT codes. It computes phonon band structures, density of states, and thermodynamic properties for crystalline materials. Widely adopted and integrated with virtually all major DFT packages, it is the de facto tool for routine phonon calculations in materials science.

**Scientific domain**: Lattice dynamics, phonons, thermodynamics, vibrational spectroscopy  
**Target user community**: Materials scientists, solid-state physicists studying vibrational properties

## Theoretical Methods
- Finite displacement method
- Harmonic approximation
- Dynamical matrix construction
- Fourier interpolation
- Born effective charges and dielectric tensors
- Non-analytical term correction (LO-TO splitting)
- Phonon group velocity
- Thermal properties (free energy, entropy, heat capacity)
- Mode Grüneisen parameters
- Quasi-harmonic approximation (QHA)
- Non-analytical term corrections (Born effective charges, dielectric tensors)
- Group theory and symmetry analysis for phonons
- Thermal properties via phonon density of states

## Capabilities (CRITICAL)
- Harmonic phonon dispersion relations
- Phonon density of states (total and projected)
- Thermal properties: free energy, entropy, heat capacity, internal energy
- Thermodynamic properties via quasi-harmonic approximation (QHA)
- Phonon group velocity
- Mode Grüneisen parameters
- Irreducible representations of phonon modes
- Band structure unfolding for supercells
- Modulation wave analysis
- Animation of phonon modes
- Interface to multiple DFT codes (VASP, Quantum ESPRESSO, ABINIT, CASTEP, CRYSTAL, DFTB+, Elk, SIESTA, TURBOMOLE, Wien2k, and more)
- Python API for custom analysis workflows
- Integration with phono3py for anharmonic extensions

**Sources**: Official phonopy documentation, GitHub repository, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - POSCAR/CONTCAR (VASP structure format) - primary
  - Various DFT code output formats via interfaces
  - FORCE_SETS (forces on displaced atoms)
  - BORN (Born effective charges and dielectric tensor)
  - mesh.yaml or band.yaml for phonon calculations
  
- **Output data types**:
  - phonon.yaml (phonon frequencies and eigenvectors)
  - thermal_properties.yaml (thermodynamic properties)
  - band.yaml (phonon band structure)
  - mesh.yaml (phonon DOS)
  - qpoints.yaml (phonon properties at specific q-points)
  - HDF5 format for large datasets
  - Animation files for visualization

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - direct integration via phonopy.structure
  - pymatgen - structure conversion utilities
  - AiiDA - workflow automation possible
  - phono3py - extends to anharmonic phonons (same author)
  
- **DFT code interfaces** (verified):
  - VASP - most common, well-tested
  - Quantum ESPRESSO - via QE output parsers
  - ABINIT - native support
  - CASTEP - .cell and .phonon formats
  - CRYSTAL - full support
  - CP2K - basic support
  - DFTB+ - supported
  - Elk - supported
  - SIESTA - supported
  - FHI-aims - supported
  - TURBOMOLE - supported
  - Wien2k - supported
  
- **Post-processing compatibility**:
  - phono3py - anharmonic phonons and thermal conductivity
  - phonopy-qha - quasi-harmonic approximation
  - phonopy-bandplot - band structure plotting
  - phononpy - Python API for custom analysis
  
- **Visualization**:
  - Built-in band structure and DOS plotting
  - Export to various formats for external tools
  - Animation generation for phonon modes

## Limitations & Known Constraints
- **Harmonic approximation**: Does not include anharmonic effects (use phono3py for anharmonicity)
- **System size**: Limited by DFT code capabilities; supercell size determines accuracy
- **Force calculation dependency**: Requires forces from DFT code; accuracy depends on DFT convergence
- **Symmetry**: Assumes perfect crystalline symmetry; defects/disorder require special treatment
- **Non-analytical terms**: Requires Born effective charges for accurate long-wavelength (Γ-point) modes in polar materials
- **Convergence**: Phonon frequencies sensitive to force constant convergence; requires careful k-point and supercell convergence tests
- **Computational cost**: Scales with supercell size; large systems require significant DFT resources

## Verification & Sources
**Primary sources**:
1. Official documentation: https://phonopy.github.io/phonopy/
2. GitHub repository: https://github.com/phonopy/phonopy
3. A. Togo and I. Tanaka, Scr. Mater. 108, 1-5 (2015) - Phonopy paper
4. Phonopy interface documentation: https://phonopy.github.io/phonopy/interfaces.html

**Secondary sources**:
1. phono3py documentation: https://phonopy.github.io/phono3py/
2. ASE phonon tutorials: Recommend phonopy for production calculations
3. Materials Project workflows: Uses phonopy for phonon calculations
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (GitHub issues, mailing list)
- Academic citations: >2,500 (Google Scholar)
- Interface support: Verified for 12+ DFT codes
