# WIEN2k-Topo (CherN/wcc modules)

## Official Resources
- Publication: Comput. Phys. Commun. (2023), doi:10.1016/j.cpc.2023.108836
- arXiv: https://arxiv.org/abs/2303.16306
- Data: Zenodo repository
- License: Check with authors

## Overview
WIEN2k-Topo comprises two Python modules (CherN.py and wcc.py) that expand the functionalities of the all-electron full-potential WIEN2k package for computing Chern and Z2 topological invariants. These complement the WloopPHI module to provide a complete toolkit for topological characterization.

**Scientific domain**: Topological invariants, all-electron calculations, band topology
**Target user community**: WIEN2k users studying topological materials

## Theoretical Methods
- Chern number calculation via Berry curvature integration
- Z2 invariant from hybrid Wannier charge centers (WCC)
- Wilson loop evolution tracking
- Fu-Kane parity method for inversion-symmetric systems
- All-electron LAPW accuracy

## Capabilities (CRITICAL)
- **Chern Number**: Berry curvature integration method
- **Z2 Invariant**: WCC evolution method
- **WIEN2k Native**: Direct interface with WIEN2k
- **All-Electron**: Full-potential LAPW accuracy
- **Validated**: Tested on known topological materials
- **WloopPHI Compatible**: Works with existing WIEN2k topology tools

**Sources**: CPC publication, arXiv preprint

## Key Strengths

### All-Electron Accuracy:
- Full-potential LAPW method
- Core electron treatment
- High accuracy for heavy elements
- Spin-orbit coupling included

### Complete Toolkit:
- Combines with WloopPHI
- Chern and Z2 together
- Multiple methods available
- Cross-validation possible

### WIEN2k Integration:
- Native WIEN2k interface
- Uses WIEN2k wavefunctions
- Consistent workflow
- Established DFT package

## Inputs & Outputs
- **Input formats**:
  - WIEN2k case files
  - Band structure data
  - k-mesh specification
  
- **Output data types**:
  - Chern numbers
  - Z2 invariants
  - WCC evolution plots
  - Berry curvature data

## Installation
```bash
# Obtain from Zenodo repository
# Place in WIEN2k SRC directory
# Configure for your WIEN2k installation
```

## Usage Examples
```bash
# After WIEN2k SCF calculation
# Calculate Chern number
python CherN.py -case material -kx 20 -ky 20

# Calculate Z2 invariant via WCC
python wcc.py -case material -plane kz0 -nk 50
```

## Performance Characteristics
- **Speed**: Depends on k-mesh density
- **Accuracy**: All-electron precision
- **Scalability**: Handles complex materials

## Limitations & Known Constraints
- **WIEN2k-specific**: Requires WIEN2k installation
- **k-mesh density**: Dense mesh needed for accuracy
- **Computational cost**: All-electron is expensive

## Comparison with Other Tools
- **vs Z2Pack**: WIEN2k-Topo all-electron, Z2Pack general
- **vs WannierTools**: Different DFT backends
- **Unique strength**: All-electron accuracy with WIEN2k

## Application Areas
- Topological insulator identification
- Chern insulator classification
- Heavy element topological materials
- Magnetic topological phases
- Surface state prediction

## Best Practices
- Use fine k-mesh for accurate invariants
- Verify with multiple methods
- Check convergence with mesh density
- Validate against known materials

## Community and Support
- CPC publication
- Zenodo data repository
- Academic correspondence

## Verification & Sources
**Primary sources**:
1. Comput. Phys. Commun. (2023), doi:10.1016/j.cpc.2023.108836
2. arXiv: https://arxiv.org/abs/2303.16306
3. Zenodo: Data and code repository

**Confidence**: VERIFIED - Published in CPC

**Verification status**: ✅ VERIFIED
- Publication: CPC peer-reviewed
- arXiv preprint: ACCESSIBLE
- Data: Zenodo repository
- Method: Chern/Z2 for WIEN2k
