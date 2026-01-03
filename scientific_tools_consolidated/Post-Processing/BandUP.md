# BandUP (Band Unfolding code for Plane-wave based calculations)

## Official Resources
- Homepage: https://bandup.readthedocs.io/
- Documentation: https://bandup.readthedocs.io/en/latest/
- Source Repository: https://github.com/band-unfolding/bandup
- License: GNU General Public License v3.0

## Overview
BandUP is a code for unfolding the electronic band structure of supercells into the primitive cell Brillouin zone. This technique is essential for analyzing calculations of impurities, defects, alloys, or surfaces, where the translational symmetry of the primitive cell is broken or obscured by the use of a supercell. BandUP recovers the effective band structure (spectral weight) in the primitive representation, making it comparable to ARPES experiments and standard band structures.

**Scientific domain**: Band unfolding, electronic structure analysis, defects, alloys  
**Target user community**: DFT users (VASP, Quantum ESPRESSO, CASTEP, ABINIT)

## Theoretical Methods
- Band Unfolding technique
- Projection of supercell wavefunctions onto primitive cell states
- Spectral weight calculation
- Effective band structure (EBS) construction
- Unfolding of spin-polarized and non-collinear bands

## Capabilities (CRITICAL)
- Unfolding bands from supercells to primitive cells
- Handling non-orthogonal supercells
- Support for spin-orbit coupling
- Calculation of spectral weights (EBS)
- Plotting tools for visualization
- Integration with major plane-wave DFT codes

**Sources**: BandUP documentation, Phys. Rev. B 89, 041407(R) (2014)

## Inputs & Outputs
- **Input formats**: WAVECAR (VASP), wavefunction files from QE/CASTEP/ABINIT
- **Output data types**: Unfolded spectral weights, band structure data files, plotting scripts

## Interfaces & Ecosystem
- **VASP**: Native support (requires WAVECAR)
- **Quantum ESPRESSO**: Supported
- **CASTEP**: Supported
- **ABINIT**: Supported
- **Python**: BandUP is a Python package

## Workflow and Usage
1. Perform DFT calculation on the supercell (generate wavefunctions).
2. Prepare input: Specify supercell matrix and primitive cell.
3. Run BandUP: `bandup unfolded`
4. Plot results: `bandup plot`

## Performance Characteristics
- Computationally efficient (post-processing)
- Memory intensive for large WAVECAR files
- Parallelization support

## Application Areas
- Doped semiconductors (impurity bands)
- Alloys and solid solutions (disorder effects)
- Interface states
- Surface reconstructions
- Vacancies and defects

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Developed by Linköping University / Paulo V. C. Medeiros

## Verification & Sources
**Primary sources**:
1. Homepage: https://bandup.readthedocs.io/
2. GitHub: https://github.com/band-unfolding/bandup
3. Publication: P. V. C. Medeiros et al., Phys. Rev. B 89, 041407(R) (2014)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Medeiros)
- Applications: Band unfolding, supercells, defects, alloys
