# ORBKIT

## Official Resources
- Homepage: https://orbkit.github.io/
- GitHub: https://github.com/orbkit/orbkit
- Documentation: https://orbkit.github.io/
- Publication: G. Hermann et al., J. Comput. Chem. 37, 1511 (2016)
- License: GNU Lesser General Public License v3.0

## Overview
ORBKIT is a parallel Python program package for post-processing wavefunction data from quantum chemical programs. It computes grid-based quantities (molecular orbitals, electron density, electrostatic potential) and non-grid-based quantities (Mulliken charges, bond orders) from various quantum chemistry output formats.

**Scientific domain**: Wavefunction analysis, electron density, molecular orbitals
**Target user community**: Quantum chemists analyzing molecular electronic structure

## Theoretical Methods
- Molecular orbital visualization
- Electron density computation
- Electrostatic potential mapping
- Mulliken population analysis
- Bond order calculation
- Transition density analysis

## Capabilities (CRITICAL)
- **Molecular Orbitals**: Grid-based visualization
- **Electron Density**: Total and spin density
- **Population Analysis**: Mulliken charges
- **Bond Orders**: Mayer and other indices
- **Multi-Code Support**: Gaussian, ORCA, MOLPRO, GAMESS, etc.
- **Parallel Processing**: MPI and multiprocessing
- **Cube Files**: Standard visualization output

**Sources**: ORBKIT documentation, JCC publication

## Key Strengths

### Multi-Code Support:
- Gaussian, ORCA, MOLPRO
- GAMESS, Turbomole, PSI4
- Molden format support
- Consistent interface

### Parallel Processing:
- MPI parallelization
- Multiprocessing support
- Efficient grid evaluation
- Large system capable

### Python Native:
- NumPy/SciPy based
- Scriptable workflows
- Jupyter compatible
- Easy installation

## Inputs & Outputs
- **Input formats**:
  - Gaussian log/fchk files
  - ORCA output
  - Molden format
  - MOLPRO, GAMESS, Turbomole
  
- **Output data types**:
  - Cube files
  - HDF5 data
  - Mulliken populations
  - Bond orders

## Installation
```bash
pip install orbkit
# Or from source
git clone https://github.com/orbkit/orbkit.git
cd orbkit
pip install -e .
```

## Usage Examples
```python
from orbkit import read, grid, core

# Read wavefunction
qc = read.main_read('molecule.molden')

# Set up grid
grid.adjust_to_geo(qc, extend=5.0)
grid.grid_init()

# Compute electron density
rho = core.rho_compute(qc)

# Save as cube file
from orbkit import output
output.cube_creator(rho, 'density.cube', qc)
```

## Performance Characteristics
- **Speed**: Parallel grid evaluation
- **Memory**: Efficient for large grids
- **Scalability**: MPI for clusters

## Limitations & Known Constraints
- **Molecular focus**: Primarily for molecules
- **Periodic systems**: Limited support
- **File formats**: Some codes not supported
- **Documentation**: Could be more extensive

## Comparison with Other Tools
- **vs Multiwfn**: ORBKIT Python-native, Multiwfn GUI
- **vs ChemTools**: Different analysis focus
- **vs Critic2**: ORBKIT molecules, Critic2 periodic
- **Unique strength**: Python ecosystem, parallel processing

## Application Areas
- Molecular orbital visualization
- Charge distribution analysis
- Reaction mechanism studies
- Excited state analysis
- Transition density mapping

## Best Practices
- Verify wavefunction file format
- Use appropriate grid resolution
- Validate against known systems
- Leverage parallel processing

## Community and Support
- GitHub repository
- JCC publication
- Active development
- LGPL licensed

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/orbkit/orbkit
2. G. Hermann et al., J. Comput. Chem. 37, 1511 (2016)
3. Documentation: https://orbkit.github.io/

**Confidence**: VERIFIED - Published in JCC

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (LGPL-3.0)
- Academic citations: Well-cited
- Active development: Maintained
