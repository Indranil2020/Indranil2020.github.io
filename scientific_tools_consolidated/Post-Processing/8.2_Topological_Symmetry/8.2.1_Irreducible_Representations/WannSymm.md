# WannSymm

## Official Resources
- Homepage: https://github.com/ccao/WannSymm
- GitHub: https://github.com/ccao/WannSymm
- Publication: G.-X. Zhi et al., Comput. Phys. Commun. 271, 108196 (2022)
- License: Check repository

## Overview
WannSymm is a symmetry analysis code for Wannier orbitals. It analyzes the symmetry properties of Wannier functions obtained from Wannier90 calculations and can symmetrize tight-binding Hamiltonians to restore exact crystallographic symmetry that may be broken due to numerical errors.

**Scientific domain**: Wannier functions, symmetry analysis, tight-binding models
**Target user community**: Researchers using Wannier90 for tight-binding model construction

## Theoretical Methods
- Symmetry operation analysis on Wannier orbitals
- Point group representation theory
- Tight-binding Hamiltonian symmetrization
- Character table analysis
- Irreducible representation assignment

## Capabilities (CRITICAL)
- **Symmetry Analysis**: Analyze symmetry of Wannier orbitals
- **Hamiltonian Symmetrization**: Restore exact crystallographic symmetry
- **DFT Interface**: Works with VASP, WIEN2k, Quantum ESPRESSO
- **SOC Support**: Handles spin-orbit coupling calculations
- **Representation Assignment**: Assign irreps to Wannier functions
- **Wannier90 Compatible**: Direct interface with Wannier90 output

**Sources**: WannSymm documentation, CPC publication

## Key Strengths

### Multi-code Support:
- VASP interface
- WIEN2k interface
- Quantum ESPRESSO interface
- Consistent workflow

### Symmetrization:
- Fix numerical symmetry breaking
- Improve TB model quality
- Essential for topological analysis
- Better band interpolation

### Published Method:
- Peer-reviewed in CPC
- Validated methodology
- Documented algorithms
- Academic support

## Inputs & Outputs
- **Input formats**:
  - Wannier90 output files (wannier90_hr.dat)
  - DFT symmetry information
  - Crystal structure
  
- **Output data types**:
  - Symmetrized Hamiltonian
  - Symmetry analysis report
  - Irreducible representations
  - Character tables

## Installation
```bash
git clone https://github.com/ccao/WannSymm.git
cd WannSymm
# Follow compilation instructions in README
```

## Usage Examples
```bash
# Create input file WannSymm.in
# Run WannSymm
./WannSymm.x < WannSymm.in

# Input file example:
DFTcode = 'VASP'
Spinors = F
prefix = 'wannier90'
OutDir = './'
```

## Performance Characteristics
- **Speed**: Fast symmetry analysis
- **Accuracy**: Restores exact symmetry
- **Reliability**: Validated against known systems

## Limitations & Known Constraints
- **Wannier90 required**: Needs Wannier90 calculation first
- **Symmetry input**: Requires proper symmetry specification
- **Learning curve**: Understanding symmetry concepts needed

## Comparison with Other Tools
- **vs symWannier**: Different symmetrization approaches
- **vs SymClosestWannier**: WannSymm post-processing focused
- **Unique strength**: Symmetry analysis and symmetrization combined

## Application Areas
- Topological material studies
- Tight-binding model improvement
- Transport calculations
- Band interpolation
- Symmetry-constrained modeling

## Best Practices
- Verify initial Wannier functions quality
- Check symmetry operations match crystal
- Validate symmetrized bands against DFT
- Use for topological invariant calculations

## Community and Support
- GitHub repository
- CPC publication
- Academic correspondence

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ccao/WannSymm
2. G.-X. Zhi et al., Comput. Phys. Commun. 271, 108196 (2022)

**Confidence**: VERIFIED - Published in CPC

**Verification status**: ✅ VERIFIED
- Official repository: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (GitHub)
- Academic citations: CPC publication
