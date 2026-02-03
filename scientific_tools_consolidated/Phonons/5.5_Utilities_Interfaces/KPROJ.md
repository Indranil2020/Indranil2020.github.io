# KPROJ

## Official Resources
- Homepage: https://arxiv.org/abs/2410.10910
- Documentation: arXiv:2410.10910
- License: Open Source

## Overview
KPROJ is a program for unfolding electronic and phononic band structures of materials modeled by supercells. It enables visualization of effective primitive cell band structures from supercell calculations, essential for studying defects, alloys, and interfaces.

**Scientific domain**: Band unfolding, supercell calculations, phonon dispersions  
**Target user community**: Researchers working with supercell calculations

## Theoretical Methods
- Band structure unfolding
- Spectral weight calculation
- Supercell to primitive cell mapping
- k-point projection
- Phonon unfolding
- Electronic unfolding

## Capabilities (CRITICAL)
- Electronic band unfolding
- Phonon band unfolding
- Supercell calculations
- Layer projection
- Twisted bilayer support
- DFT code interfaces
- Spectral function output

## Key Strengths

### Dual Unfolding:
- Both electronic and phononic
- Unified framework
- Consistent methodology
- Versatile application

### Supercell Support:
- Defect calculations
- Alloy studies
- Interface systems
- Twisted structures

## Inputs & Outputs
- **Input formats**:
  - DFT output files
  - Phonon eigenvectors
  - Supercell information
  
- **Output data types**:
  - Unfolded band structures
  - Spectral weights
  - Effective dispersions

## Interfaces & Ecosystem
- **VASP**: Primary interface
- **Phonopy**: Phonon input
- **Python**: Analysis tools


## Advanced Features
- **Dual unfolding**: Both electronic and phononic band structures
- **Layer projection**: Layer-resolved spectral weights
- **Twisted bilayer support**: Moiré superlattice calculations
- **Spectral weight calculation**: Bloch character analysis
- **Multiple DFT codes**: VASP and Phonopy interfaces

## Performance Characteristics
- Post-processing tool: Fast
- Scales with supercell size
- Python-based implementation

## Computational Cost
- Supercell DFT/phonon calculation: External (dominant cost)
- KPROJ unfolding: Fast (minutes)
- Overall: Minimal overhead for unfolding

## Best Practices
- Use commensurate supercell-primitive mapping
- Check spectral weight convergence
- Validate against primitive cell bands when possible
- Use appropriate broadening for visualization

## Limitations & Known Constraints
- Recent development
- Limited documentation
- Specific DFT codes
- Expert-level tool

## Application Areas
- Defect phonons
- Alloy band structures
- Interface studies
- Twisted bilayers
- Disordered systems

## Verification & Sources
**Primary sources**:
1. arXiv: https://arxiv.org/abs/2410.10910

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Published methodology
- arXiv preprint
