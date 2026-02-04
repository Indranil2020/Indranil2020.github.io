# ChemTools

## Official Resources
- Homepage: https://chemtools.org/
- GitHub: https://github.com/theochem/chemtools
- Documentation: https://chemtools.org/
- PyPI: https://pypi.org/project/chemtools/
- License: GNU General Public License v3.0

## Overview
ChemTools is a Python library for interpreting the results of quantum chemistry calculations using conceptual density functional theory (CDFT). It provides tools for computing and visualizing chemical reactivity descriptors including Fukui functions, dual descriptor, electrophilicity, and other local and global reactivity indicators.

**Scientific domain**: Conceptual DFT, chemical reactivity, Fukui functions
**Target user community**: Computational chemists studying reactivity and selectivity

## Theoretical Methods
- Conceptual density functional theory
- Fukui functions (f+, f-, f0)
- Dual descriptor
- Global reactivity descriptors (chemical potential, hardness)
- Local softness and electrophilicity
- Condensed-to-atoms descriptors

## Capabilities (CRITICAL)
- **Fukui Functions**: Nucleophilic, electrophilic, radical
- **Dual Descriptor**: Electrophilic/nucleophilic regions
- **Global Descriptors**: Hardness, softness, electronegativity
- **Local Descriptors**: Local softness, electrophilicity index
- **Visualization**: 3D isosurface generation
- **Multiple Formats**: Gaussian, ORCA, wfn/wfx support

**Sources**: ChemTools documentation, theochem group

## Key Strengths

### Conceptual DFT Focus:
- Comprehensive CDFT implementation
- All major descriptors included
- Theoretical rigor
- Active development

### Python Native:
- pip installable
- NumPy/SciPy based
- Jupyter notebook support
- Scriptable workflows

### Visualization:
- 3D isosurface generation
- VMD integration
- Publication-quality output
- Interactive analysis

## Inputs & Outputs
- **Input formats**:
  - Gaussian fchk files
  - wfn/wfx wavefunction files
  - ORCA output
  - Molden files
  
- **Output data types**:
  - Reactivity descriptors
  - Cube files for visualization
  - Condensed values
  - Isosurfaces

## Installation
```bash
pip install chemtools
# Or from source
git clone https://github.com/theochem/chemtools.git
cd chemtools
pip install -e .
```

## Usage Examples
```python
from chemtools import LocalConceptualDFT, UniformGrid

# Load molecule calculations (neutral, cation, anion)
mol = LocalConceptualDFT.from_file(['neutral.fchk', 'cation.fchk', 'anion.fchk'])

# Compute Fukui functions on a grid
grid = UniformGrid.from_molecule(mol.molecule)
fukui_plus = mol.fukui_function(grid.points, 'plus')
fukui_minus = mol.fukui_function(grid.points, 'minus')
dual = mol.dual_descriptor(grid.points)

# Generate cube file
grid.generate_cube('dual_descriptor.cube', dual)
```

## Performance Characteristics
- **Speed**: Efficient grid-based evaluation
- **Memory**: Depends on grid resolution
- **Accuracy**: Quantum chemistry precision

## Limitations & Known Constraints
- **Three calculations**: Requires N, N+1, N-1 electron calculations
- **File formats**: Limited DFT code support
- **Learning curve**: CDFT concepts required
- **Documentation**: Evolving

## Comparison with Other Tools
- **vs Multiwfn**: ChemTools Python-native, Multiwfn GUI
- **vs ADF reactivity**: ChemTools open-source
- **vs Gaussian population**: ChemTools more comprehensive
- **Unique strength**: Dedicated CDFT library, Python ecosystem

## Application Areas
- Regioselectivity prediction
- Reaction mechanism analysis
- Electrophilic/nucleophilic site identification
- Catalyst design
- Drug-receptor interactions

## Best Practices
- Use consistent basis sets for N, N±1
- Verify SCF convergence for all states
- Compare with chemical intuition
- Validate on known reactions

## Community and Support
- GitHub repository
- theochem research group
- Active development
- GPL v3 licensed

## Verification & Sources
**Primary sources**:
1. Homepage: https://chemtools.org/
2. GitHub: https://github.com/theochem/chemtools
3. P. W. Ayers et al. (theochem group publications)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, GPL-3.0)
- Developer: theochem group
- Active development: Maintained
