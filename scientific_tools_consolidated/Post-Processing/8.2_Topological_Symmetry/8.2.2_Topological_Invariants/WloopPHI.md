# WloopPHI

## Official Resources
- Publication: H. Saini et al., Comput. Phys. Commun. 270, 108147 (2022)
- arXiv: https://arxiv.org/abs/2008.08124
- License: Check with authors

## Overview
WloopPHI is a Python code that expands the features of WIEN2k by enabling characterization of Weyl semimetals through Wilson loop calculations. It computes the winding of hybrid Wannier charge centers to identify and characterize Weyl points, nodal lines, and other topological features in the band structure.

**Scientific domain**: Weyl semimetals, topological semimetals, Wilson loop analysis
**Target user community**: WIEN2k users studying topological semimetals and nodal structures

## Theoretical Methods
- Wilson loop calculations
- Hybrid Wannier charge centers (WCC)
- Winding number computation
- Berry phase integration
- Weyl point chirality determination
- Nodal line characterization

## Capabilities (CRITICAL)
- **WIEN2k Integration**: Direct interface with WIEN2k
- **Weyl Points**: Chirality and position determination
- **Nodal Lines**: Nodal structure characterization
- **Wilson Loops**: WCC evolution calculation
- **Chern Numbers**: Surface Chern number computation
- **Berry Phase**: Integrated Berry phase calculation

**Sources**: CPC publication, arXiv preprint

## Key Strengths

### WIEN2k Native:
- Direct WIEN2k interface
- Uses WIEN2k wavefunctions
- Consistent with LAPW method
- All-electron accuracy

### Semimetal Focus:
- Weyl point identification
- Chirality determination
- Nodal line detection
- Fermi arc prediction

### Validated Implementation:
- Tested on TaAs, Na3Bi, CaAgAs
- Published benchmarks
- Verified methodology
- Robust algorithms

## Inputs & Outputs
- **Input formats**:
  - WIEN2k case files
  - Band structure data
  - Wannier functions
  
- **Output data types**:
  - Wilson loop spectra
  - Winding numbers
  - Weyl point positions
  - Chirality values

## Installation
```bash
# Obtain from authors or CPC library
# Requires WIEN2k installation
python setup.py install
```

## Usage Examples
```bash
# After WIEN2k SCF calculation
# Generate Wilson loop on a sphere around suspected Weyl point
python wloopphi.py -case material -kpath sphere -center 0.5 0.0 0.0 -radius 0.1
```

## Performance Characteristics
- **Speed**: Depends on k-point density
- **Accuracy**: All-electron LAPW precision
- **Scalability**: Handles complex band structures

## Limitations & Known Constraints
- **WIEN2k-specific**: Only works with WIEN2k
- **Manual setup**: Requires careful k-path definition
- **Computational cost**: Dense k-mesh needed for accuracy

## Comparison with Other Tools
- **vs WannierTools**: WloopPHI WIEN2k-based, WannierTools Wannier90
- **vs Z2Pack**: Different DFT backends
- **Unique strength**: Native WIEN2k integration, all-electron accuracy

## Application Areas
- Weyl semimetal identification
- Dirac semimetal characterization
- Nodal line semimetals
- Topological metal analysis
- Fermi arc prediction

## Best Practices
- Use fine k-mesh around suspected nodes
- Verify Weyl point positions with band structure
- Check chirality consistency
- Validate with known materials

## Community and Support
- CPC publication
- Academic correspondence
- WIEN2k community

## Verification & Sources
**Primary sources**:
1. H. Saini et al., Comput. Phys. Commun. 270, 108147 (2022)
2. arXiv: https://arxiv.org/abs/2008.08124

**Confidence**: VERIFIED - Published in CPC

**Verification status**: ✅ VERIFIED
- Publication: CPC peer-reviewed
- arXiv preprint: ACCESSIBLE
- Method: Wilson loop for WIEN2k
- Benchmarks: TaAs, Na3Bi, CaAgAs validated
