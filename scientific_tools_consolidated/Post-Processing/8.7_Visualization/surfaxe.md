# surfaxe

## Official Resources
- Source Repository: https://github.com/SMTG-Bham/surfaxe
- Documentation: https://surfaxe.readthedocs.io/
- PyPI: https://pypi.org/project/surfaxe/
- License: MIT License

## Overview
**surfaxe** is a Python package for dealing with slabs for first principles calculations. It automates surface energy convergence checks, work function calculations, bond length analysis, and slab thickness convergence for surface science DFT calculations with VASP.

**Scientific domain**: Surface science, slab calculations, surface energy, work function  
**Target user community**: Researchers performing DFT surface calculations and needing automated convergence and analysis

## Theoretical Methods
- Surface energy calculation
- Work function from electrostatic potential
- Slab convergence analysis
- Bond length analysis vs depth
- Surface reconstruction analysis
- VASP output parsing

## Capabilities (CRITICAL)
- Surface energy calculation and convergence
- Work function calculation
- Average bond length vs slab depth
- Surface energy convergence with slab thickness
- Electrostatic potential analysis
- VASP output parsing (LOCPOT, OUTCAR)
- Automatic slab generation helpers
- Plotting and visualization

**Sources**: GitHub repository, JOSS

## Key Strengths

### Automated Surface Analysis:
- Surface energy convergence
- Work function from LOCPOT
- Bond length depth profiles
- No manual data extraction

### VASP Integration:
- LOCPOT parsing for potential
- OUTCAR parsing for energies
- POSCAR structure handling
- Consistent with VASP workflow

### Convergence:
- Slab thickness convergence
- k-point convergence for surfaces
- Vacuum size convergence
- Systematic testing

## Inputs & Outputs
- **Input formats**:
  - VASP output files (LOCPOT, OUTCAR, POSCAR)
  - Surface specification
  
- **Output data types**:
  - Surface energies
  - Work functions
  - Bond length profiles
  - Convergence plots

## Interfaces & Ecosystem
- **VASP**: Primary DFT backend
- **pymatgen**: Structure handling
- **Matplotlib**: Visualization
- **Python**: Scripting

## Performance Characteristics
- **Speed**: Fast (post-processing)
- **Accuracy**: DFT-level
- **System size**: Any slab size
- **Automation**: High

## Computational Cost
- **Analysis**: Seconds
- **VASP pre-requisite**: Hours (separate)
- **Typical**: Very efficient

## Limitations & Known Constraints
- **VASP only**: No QE or other code support
- **Surface only**: No bulk analysis
- **3D slabs**: No 2D material specific features
- **Limited to specific analyses**: Surface-focused

## Comparison with Other Codes
- **vs VESTA**: surfaxe calculates surface properties, VESTA visualizes
- **vs pymatgen**: surfaxe is surface-specialized, pymatgen is general
- **vs ASE**: surfaxe has surface-specific analysis, ASE is general
- **Unique strength**: Automated surface energy convergence, work function, and slab analysis for VASP

## Application Areas

### Surface Science:
- Surface energy determination
- Work function calculation
- Surface stability analysis
- Surface reconstruction

### Catalysis:
- Catalyst surface properties
- Adsorption site preparation
- Surface energy for stability
- Work function for charge transfer

### Semiconductor Surfaces:
- Fermi level pinning
- Surface band bending
- Interface properties
- Schottky barrier estimation

### 2D Materials:
- Exfoliation energy
- Surface energy of layers
- Work function of monolayers
- Van der Waals surfaces

## Best Practices

### Slab Convergence:
- Test convergence with slab thickness
- Use sufficient vacuum
- Check k-point convergence
- Monitor surface energy convergence

### Work Function:
- Use well-converged LOCPOT
- Average over vacuum region
- Check dipole correction
- Compare with experiment

## Community and Support
- Open source (MIT License)
- PyPI installation available
- ReadTheDocs documentation
- Developed at University of Birmingham (SMTG)
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/SMTG-Bham/surfaxe
2. Documentation: https://surfaxe.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- PyPI: AVAILABLE
- Active development: Ongoing
- Specialized strength: Automated surface energy convergence, work function, and slab analysis for VASP
