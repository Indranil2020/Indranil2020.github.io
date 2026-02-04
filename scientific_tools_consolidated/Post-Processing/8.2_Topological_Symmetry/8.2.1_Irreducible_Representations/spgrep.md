# spgrep

## Official Resources
- Homepage: https://spglib.github.io/spgrep/
- GitHub: https://github.com/spglib/spgrep
- Documentation: https://spglib.github.io/spgrep/
- PyPI: https://pypi.org/project/spgrep/
- Publication: K. Shinohara et al., J. Open Source Softw. 8, 5269 (2023)
- License: BSD 3-Clause License

## Overview
spgrep is a Python package for on-the-fly generation of space-group irreducible representations. It computes irreducible representations (irreps) and their characters for any space group at arbitrary k-points, without relying on pre-tabulated databases. This enables flexible symmetry analysis for band structure calculations.

**Scientific domain**: Crystallographic symmetry, irreducible representations, electronic structure
**Target user community**: Computational materials scientists, condensed matter physicists

## Theoretical Methods
- Induction from point group representations
- Projective representations for non-symmorphic groups
- On-the-fly irrep construction
- Compatibility with spglib symmetry operations
- Spin representations (double groups)

## Capabilities (CRITICAL)
- **On-the-fly Irrep Generation**: No pre-tabulated database needed
- **Arbitrary k-points**: Works at any k-vector
- **Double Groups**: Spin-orbit coupling support
- **Magnetic Groups**: Magnetic space group support
- **Character Tables**: Automatic generation
- **Representation Matrices**: Full matrix representations
- **spglib Integration**: Uses spglib for symmetry operations

**Sources**: spgrep documentation, JOSS publication

## Key Strengths

### On-the-fly Computation:
- No lookup tables required
- Works for any k-point
- Flexible and general
- Always up-to-date with spglib

### Python Native:
- Easy installation via pip
- NumPy-based computations
- Jupyter notebook friendly
- Scriptable workflows

### Modern Implementation:
- Clean API design
- Well-documented
- JOSS peer-reviewed
- Active maintenance

## Inputs & Outputs
- **Input formats**:
  - spglib symmetry operations
  - k-point coordinates
  - Crystal structure
  
- **Output data types**:
  - Irreducible representations
  - Character tables
  - Representation matrices
  - Irrep labels

## Installation
```bash
pip install spgrep
```

## Usage Examples
```python
import spgrep
from spgrep import get_spacegroup_irreps
import spglib

# Get symmetry operations
cell = (lattice, positions, numbers)
symmetry = spglib.get_symmetry(cell)

# Generate irreps at Gamma point
kpoint = [0, 0, 0]
irreps = get_spacegroup_irreps(symmetry, kpoint)

# Get character table
for irrep in irreps:
    print(f"Irrep: {irrep.label}, Dimension: {irrep.dimension}")
```

## Performance Characteristics
- **Speed**: Fast on-the-fly generation
- **Memory**: Efficient, no large databases
- **Flexibility**: Works for any k-point

## Limitations & Known Constraints
- **spglib dependency**: Requires spglib for symmetry
- **Numerical precision**: Floating-point k-points
- **Convention**: Uses spglib conventions

## Comparison with Other Tools
- **vs SpaceGroupIrep**: spgrep Python, SpaceGroupIrep Mathematica
- **vs IrRep**: Different approaches to irrep computation
- **vs phonopy**: spgrep general purpose, phonopy phonon-focused
- **Unique strength**: On-the-fly computation without databases

## Application Areas
- Band structure symmetry analysis
- Phonon symmetry classification
- Selection rules derivation
- Topological band theory
- Crystal field splitting

## Best Practices
- Verify k-point conventions match your DFT code
- Use with spglib for consistent symmetry
- Check irrep dimensions for validation
- Consider numerical tolerance for k-points

## Community and Support
- GitHub issue tracker
- JOSS publication for citation
- Part of spglib ecosystem
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/spglib/spgrep
2. Documentation: https://spglib.github.io/spgrep/
3. K. Shinohara et al., J. Open Source Softw. 8, 5269 (2023)

**Confidence**: VERIFIED - Published in JOSS

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, BSD-3)
- Developer: spglib team
- Academic citations: JOSS publication
- Active development: Regular releases
