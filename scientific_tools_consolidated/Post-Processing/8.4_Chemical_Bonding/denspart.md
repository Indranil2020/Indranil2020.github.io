# denspart

## Official Resources
- Homepage: https://github.com/theochem/denspart
- GitHub: https://github.com/theochem/denspart
- Documentation: GitHub README
- License: GNU General Public License v3.0

## Overview
denspart is a Python package for atoms-in-molecules density partitioning. It implements various stockholder-type partitioning schemes including Iterative Stockholder Atoms (ISA), Minimal Basis Iterative Stockholder (MBIS), and Gaussian Iterative Stockholder (GISA) methods for computing atomic charges and multipole moments.

**Scientific domain**: Charge partitioning, atoms-in-molecules, stockholder methods
**Target user community**: Computational chemists developing force fields and analyzing charge distributions

## Theoretical Methods
- Iterative Stockholder Atoms (ISA)
- Minimal Basis Iterative Stockholder (MBIS)
- Gaussian Iterative Stockholder (GISA)
- Multipole moment calculation
- Stockholder partitioning

## Capabilities (CRITICAL)
- **ISA Method**: Iterative stockholder partitioning
- **MBIS Method**: Minimal basis stockholder
- **GISA Method**: Gaussian basis stockholder
- **Multipole Moments**: Full multipole expansion
- **Python Native**: Modern Python implementation
- **NumPy Based**: Efficient array operations

**Sources**: denspart GitHub repository, theochem group

## Key Strengths

### Modern Implementation:
- Python 3 native
- NumPy/SciPy based
- Clean API design
- Active development

### Multiple Methods:
- ISA, MBIS, GISA
- Consistent interface
- Method comparison
- Flexibility

### theochem Ecosystem:
- IOData integration
- Grid library support
- ChemTools compatible
- Consistent tools

## Inputs & Outputs
- **Input formats**:
  - NumPy .npz density files
  - Grid-based density data
  - Via IOData parsers
  
- **Output data types**:
  - Atomic charges
  - Multipole moments
  - Partitioned densities
  - NumPy arrays

## Installation
```bash
pip install denspart
# Or from source
git clone https://github.com/theochem/denspart.git
cd denspart
pip install -e .
```

## Usage Examples
```python
from denspart.mbis import MBISProModel
import numpy as np

# Load density data
data = np.load('density.npz')

# Perform MBIS partitioning
pro_model = MBISProModel.from_geometry(atnums, atcoords)
pro_model.run(density, grid)

# Get atomic charges
charges = pro_model.charges
multipoles = pro_model.multipole_moments
```

## Performance Characteristics
- **Speed**: Efficient NumPy operations
- **Memory**: Scales with grid size
- **Convergence**: Iterative methods

## Limitations & Known Constraints
- **Grid input**: Requires pre-computed density grid
- **Documentation**: Evolving
- **Molecular focus**: Primarily for molecules
- **Dependencies**: Requires theochem ecosystem

## Comparison with Other Tools
- **vs DDEC**: denspart open-source, different methods
- **vs Hirshfeld**: denspart implements iterative variants
- **vs horton-part**: Related but separate packages
- **Unique strength**: Modern Python implementation of ISA/MBIS

## Application Areas
- Force field parameterization
- Charge distribution analysis
- Electrostatic property prediction
- Reactivity indices
- Machine learning features

## Best Practices
- Use appropriate grid density
- Compare multiple methods
- Verify charge conservation
- Validate against known systems

## Community and Support
- GitHub repository
- theochem research group
- GPL v3 licensed
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/theochem/denspart
2. theochem group (McMaster University)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Source code: OPEN (GPL-3.0)
- Developer: theochem group
- Method: ISA/MBIS implementation
