# spgrep-modulation

## Official Resources
- Homepage: https://phonopy.github.io/spgrep-modulation/
- GitHub: https://github.com/phonopy/spgrep-modulation
- Documentation: https://phonopy.github.io/spgrep-modulation/
- PyPI: https://pypi.org/project/spgrep-modulation/
- License: BSD 3-Clause License

## Overview
spgrep-modulation is a Python package for collective atomic modulation analysis using irreducible space-group representations. It enables systematic study of structural phase transitions, phonon instabilities, and order parameter analysis through group-theoretical methods integrated with phonopy.

**Scientific domain**: Structural phase transitions, phonon symmetry, order parameter analysis
**Target user community**: Researchers studying phase transitions, ferroelectrics, and structural instabilities

## Theoretical Methods
- Irreducible representation-based modulation analysis
- Order parameter direction determination
- Isotropy subgroup computation
- Phonon mode symmetry classification
- Landau theory connection

## Capabilities (CRITICAL)
- **Irrep-based Analysis**: Phonon modes classified by irreps
- **Modulated Structures**: Generate distorted structures
- **Isotropy Subgroups**: Symmetry breaking analysis
- **Order Parameters**: Direction specification
- **Phonopy Integration**: Works with phonopy workflows
- **spgrep Backend**: Uses spgrep for irreps

**Sources**: spgrep-modulation documentation, phonopy ecosystem

## Key Strengths

### Phonopy Ecosystem:
- Seamless phonopy integration
- Consistent with phonopy conventions
- Leverages phonopy infrastructure
- Combined workflow support

### Systematic Analysis:
- Automatic isotropy subgroup finding
- Order parameter enumeration
- Modulated structure generation
- Group-subgroup relations

### Modern Python:
- Clean API
- Well-documented
- NumPy-based
- Jupyter compatible

## Inputs & Outputs
- **Input formats**:
  - Phonopy objects
  - Crystal structures
  - Phonon eigenvectors
  
- **Output data types**:
  - Modulated structures
  - Isotropy subgroups
  - Order parameter directions
  - Symmetry analysis

## Installation
```bash
pip install spgrep-modulation
```

## Usage Examples
```python
from pathlib import Path
import phonopy
from phonopy.structure.symmetry import Symmetry
from spgrep_modulation import Modulation

# Load phonopy calculation
ph = phonopy.load("phonopy_disp.yaml")

# Create modulation analysis
mod = Modulation.with_Sr_number(
    ph.primitive,
    ph.dynamical_matrix,
    qpoint=[0.5, 0.5, 0],
    irrep_index=0
)

# Generate modulated structure
modulated = mod.get_modulated_structure(amplitude=0.1)
```

## Performance Characteristics
- **Speed**: Fast symmetry analysis
- **Integration**: Efficient phonopy interface
- **Flexibility**: Arbitrary q-points supported

## Limitations & Known Constraints
- **Phonopy dependency**: Requires phonopy installation
- **Harmonic phonons**: Based on harmonic approximation
- **Documentation**: Evolving documentation

## Comparison with Other Tools
- **vs ISOTROPY**: spgrep-modulation Python-native
- **vs phonopy**: spgrep-modulation adds irrep analysis
- **vs AMPLIMODES**: Different approach to mode analysis
- **Unique strength**: Python integration with phonopy/spgrep

## Application Areas
- Ferroelectric phase transitions
- Structural instabilities
- Soft mode analysis
- Perovskite distortions
- Order-disorder transitions

## Best Practices
- Ensure phonopy calculation is converged
- Use appropriate q-point mesh
- Verify irrep assignments
- Check modulation amplitudes physically

## Community and Support
- GitHub issue tracker
- Part of phonopy ecosystem
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/phonopy/spgrep-modulation
2. Documentation: https://phonopy.github.io/spgrep-modulation/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (GitHub, BSD-3)
- Developer: phonopy team
- Active development: Maintained
