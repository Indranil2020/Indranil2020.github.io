# DensityTool

## Official Resources
- Source Repository: https://github.com/llodeiro/DensityTool
- Documentation: Included in repository
- License: Open source

## Overview
**DensityTool** is a FORTRAN post-processing program for space- and spin-resolved density of states from VASP. It provides detailed spatial and spin decomposition of the DOS, enabling atom-, shell-, and orbital-projected DOS analysis beyond standard VASP output.

**Scientific domain**: DOS analysis, spin-resolved density of states, VASP post-processing  
**Target user community**: Researchers needing detailed spatial and spin decomposition of DOS from VASP calculations

## Theoretical Methods
- Space-resolved density of states
- Spin-resolved DOS analysis
- Atom-projected DOS
- Shell-projected DOS (s, p, d, f)
- Orbital-projected DOS
- VASP PROCAR/DOSCAR parsing

## Capabilities (CRITICAL)
- Space-resolved DOS calculation
- Spin-resolved DOS decomposition
- Atom-projected DOS
- Shell (s, p, d, f) projected DOS
- Orbital-projected DOS
- VASP output parsing
- Custom region DOS integration

**Sources**: GitHub repository, Comput. Phys. Commun. 277, 108384 (2022)

## Key Strengths

### Space-Resolved DOS:
- DOS in custom spatial regions
- Not limited to atomic spheres
- Integration over user-defined volumes
- Unique spatial decomposition

### Spin Resolution:
- Full spin-resolved analysis
- Up/down spin channels
- Spin density decomposition
- Magnetic DOS analysis

### Beyond Standard VASP:
- More detailed than standard PROCAR
- Custom spatial regions
- Shell and orbital decomposition
- Flexible output

## Inputs & Outputs
- **Input formats**:
  - VASP PROCAR
  - VASP DOSCAR
  - Spatial region definitions
  
- **Output data types**:
  - Space-resolved DOS
  - Spin-resolved DOS
  - Atom/shell/orbital projected DOS
  - Custom region DOS

## Interfaces & Ecosystem
- **VASP**: Primary DFT backend
- **FORTRAN**: Core computation

## Performance Characteristics
- **Speed**: Fast (post-processing)
- **Accuracy**: VASP-level
- **System size**: Limited by PROCAR size
- **Memory**: Moderate

## Computational Cost
- **DOS analysis**: Seconds to minutes
- **VASP pre-requisite**: Hours (separate)
- **Typical**: Very efficient

## Limitations & Known Constraints
- **VASP only**: No QE or other code support
- **FORTRAN**: Compilation required
- **PROCAR dependency**: Requires LORBIT=11
- **Limited documentation**: Could be more extensive

## Comparison with Other Codes
- **vs pyprocar**: DensityTool has space-resolved DOS, pyprocar is more general
- **vs VASPKIT**: DensityTool has spatial decomposition, VASPKIT is broader
- **vs sumo**: DensityTool has unique space-resolved DOS
- **Unique strength**: Space-resolved and spin-resolved DOS decomposition from VASP, custom spatial regions

## Application Areas

### Magnetic Materials:
- Spin-resolved DOS
- Magnetic moment analysis
- Exchange splitting
- Spin-polarized surface states

### Surface Science:
- Surface-localized DOS
- Interface state analysis
- Vacuum region DOS
- Spatial DOS mapping

### Defect Analysis:
- Defect state localization
- Spatial extent of defect states
- Spin polarization at defects
- Impurity state decomposition

## Best Practices

### VASP Setup:
- Set LORBIT=11 for PROCAR
- Use sufficient k-points for DOS
- Include spin polarization if needed
- Use appropriate RWIGS

### Spatial Analysis:
- Define meaningful spatial regions
- Test convergence with grid
- Compare with atomic sphere results
- Validate against total DOS

## Community and Support
- Open source on GitHub
- Published in Comput. Phys. Commun.
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/llodeiro/DensityTool
2. L. Lodeiro and T. Rauch, Comput. Phys. Commun. 277, 108384 (2022)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Published methodology: Comput. Phys. Commun.
- Specialized strength: Space-resolved and spin-resolved DOS decomposition from VASP, custom spatial regions
