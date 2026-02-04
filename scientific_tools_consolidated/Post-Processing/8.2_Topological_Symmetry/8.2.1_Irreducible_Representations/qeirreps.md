# qeirreps

## Official Resources
- Homepage: https://github.com/mizoguche/qeirreps
- GitHub: https://github.com/mizoguche/qeirreps
- Publication: A. Matsugatani et al., Comput. Phys. Commun. 264, 107948 (2021)
- License: MIT License

## Overview
qeirreps is an open-source program that computes irreducible representations of Bloch wavefunctions from Quantum ESPRESSO output. It analyzes the symmetry character of electronic bands at high-symmetry k-points, enabling topological classification and symmetry-based materials analysis.

**Scientific domain**: Electronic structure symmetry, irreducible representations, topological materials
**Target user community**: Quantum ESPRESSO users studying band topology and symmetry properties

## Theoretical Methods
- Bloch wavefunction symmetry analysis
- Irreducible representation assignment
- Non-symmorphic space group handling
- Spinor representations (spin-orbit coupling)
- Symmetry eigenvalue computation

## Capabilities (CRITICAL)
- **QE Integration**: Direct interface with Quantum ESPRESSO
- **Irrep Assignment**: Automatic band symmetry classification
- **Spin-Orbit**: Handles SOC calculations
- **Non-symmorphic**: Proper phase factors for glide/screw
- **Topological Analysis**: Symmetry indicator extraction
- **Multiple k-points**: Batch processing capability

**Sources**: qeirreps documentation, CPC publication

## Key Strengths

### Quantum ESPRESSO Native:
- Direct QE output parsing
- Consistent with QE conventions
- Handles QE-specific formats
- No intermediate conversion

### Topological Materials:
- Symmetry indicator computation
- Band topology diagnosis
- Compatibility with TQC database
- Fragile topology detection

### Robust Implementation:
- Validated against known materials
- Handles complex structures
- Proper SOC treatment
- Non-symmorphic accuracy

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO output files
  - Wavefunction files (wfc*.dat)
  - Band structure output
  
- **Output data types**:
  - Irreducible representation labels
  - Symmetry eigenvalues
  - Character tables
  - Topological indicators

## Installation
```bash
git clone https://github.com/mizoguche/qeirreps.git
cd qeirreps
make
```

## Usage Examples
```bash
# Run after QE band calculation
qeirreps.x -pp < input.nml

# Input namelist
&inputpp
  prefix = 'silicon'
  outdir = './tmp/'
  filband = 'bands.dat'
/
```

## Performance Characteristics
- **Speed**: Fast symmetry analysis
- **Memory**: Moderate for wavefunction handling
- **Scalability**: Handles multiple bands efficiently

## Limitations & Known Constraints
- **QE-specific**: Only works with Quantum ESPRESSO
- **k-point selection**: Requires high-symmetry points
- **Wavefunction access**: Needs QE wavefunction files

## Comparison with Other Tools
- **vs irvsp**: qeirreps for QE, irvsp for VASP
- **vs IrRep**: Both compute irreps, different DFT interfaces
- **vs spgrep**: qeirreps DFT-specific, spgrep general
- **Unique strength**: Native Quantum ESPRESSO integration

## Application Areas
- Topological insulator identification
- Symmetry-protected topology
- Band inversion analysis
- Weyl/Dirac semimetal classification
- Elementary band representation decomposition

## Best Practices
- Use consistent symmetry settings in QE
- Include enough bands for analysis
- Verify k-point high-symmetry assignment
- Check against known materials

## Community and Support
- GitHub repository
- CPC publication methodology
- Academic development (Watanabe group)

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mizoguche/qeirreps
2. A. Matsugatani et al., Comput. Phys. Commun. 264, 107948 (2021)

**Confidence**: VERIFIED - Published in CPC

**Verification status**: ✅ VERIFIED
- Official repository: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (GitHub, MIT)
- Academic citations: CPC publication
- Active development: Maintained
