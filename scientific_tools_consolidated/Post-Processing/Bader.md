# Bader (Henkelman Group)

## Official Resources
- Homepage: https://theory.cm.utexas.edu/henkelman/code/bader/
- Documentation: https://theory.cm.utexas.edu/henkelman/code/bader/
- Source Repository: https://github.com/henkelmanlab/bader
- License: GNU General Public License v2.0 (or v3.0)

## Overview
The Henkelman Group Bader code is a widely used software for performing Bader charge analysis. It partitions the charge density of a system into atomic volumes based on zero-flux surfaces of the electron density gradient. By integrating the charge within these Bader volumes, it assigns partial charges to atoms, providing a physically motivated way to define oxidation states and charge transfer.

**Scientific domain**: Charge analysis, topology of electron density, atoms in molecules  
**Target user community**: Computational chemists, materials scientists

## Theoretical Methods
- Quantum Theory of Atoms in Molecules (QTAIM)
- Zero-flux surfaces of electron density gradient
- Grid-based partitioning algorithm
- Integration of charge density
- Atomic volume calculation

## Capabilities (CRITICAL)
- Calculation of Bader charges
- Atomic volume determination
- Identification of bond critical points (optional features)
- Robust algorithm for grid-based data
- Handling of all-electron and pseudopotential densities (using core charge reconstruction)
- Output of atomic volumes for visualization

**Sources**: Bader documentation, Comp. Mater. Sci. 36, 354 (2006)

## Inputs & Outputs
- **Input formats**: CHGCAR (VASP), cube files (Gaussian, etc.), spin-density files
- **Output data types**: ACF.dat (charges/positions), BCF.dat (volume boundaries), AVF.dat (atomic volumes)

## Interfaces & Ecosystem
- **VASP**: Primary target (CHGCAR support)
- **Gaussian/Q-Chem**: Via cube file conversion
- **ASE**: Interface to run bader
- **Pymatgen**: Analysis class for parsing output

## Workflow and Usage
1. Perform DFT calculation (save charge density).
2. For PAW/Pseudopotentials: Generate reference core density (AECCAR0 + AECCAR2 in VASP).
3. Run Bader: `bader CHGCAR -ref CHGCAR_sum`
4. Read `ACF.dat` for partial charges.

## Performance Characteristics
- Fast analysis (seconds to minutes)
- Scales linearly with grid size
- Memory usage proportional to charge density grid

## Application Areas
- Oxidation state assignment
- Charge transfer in ionic/covalent bonds
- Surface adsorption analysis
- Electrocatalysis

## Community and Support
- Developed by Henkelman Group (UT Austin)
- Active support via forum
- Standard tool in VASP community

## Verification & Sources
**Primary sources**:
1. Homepage: https://theory.cm.utexas.edu/henkelman/code/bader/
2. Publication: G. Henkelman, A. Arnaldsson, and H. Jonsson, Comp. Mater. Sci. 36, 354 (2006)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Henkelman Group)
- Applications: Charge analysis, QTAIM, VASP integration
