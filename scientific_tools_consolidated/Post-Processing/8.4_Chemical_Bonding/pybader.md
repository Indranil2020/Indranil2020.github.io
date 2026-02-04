# pybader

## Official Resources
- Homepage: https://github.com/adam-kerrigan/pybader
- GitHub: https://github.com/adam-kerrigan/pybader
- PyPI: https://pypi.org/project/pybader/
- Documentation: GitHub README and examples
- License: MIT License

## Overview
pybader is a threaded Python implementation of grid-based Bader charge analysis. It provides a fast, pure-Python alternative to the original Bader code, with support for multiple file formats and integration with Python workflows for high-throughput analysis.

**Scientific domain**: Bader charge analysis, QTAIM, charge partitioning
**Target user community**: Python-based materials science workflows, high-throughput screening

## Theoretical Methods
- Bader's Quantum Theory of Atoms in Molecules (QTAIM)
- Zero-flux surface partitioning
- Near-grid and weight methods
- Charge density integration
- Basin volume calculation

## Capabilities (CRITICAL)
- **Python Native**: Pure Python implementation
- **Multi-threaded**: Parallel processing support
- **Multiple Formats**: VASP, CHGCAR, cube files
- **Bader Volumes**: Atomic basin volumes
- **Charge Integration**: Accurate charge partitioning
- **Refinement Options**: Multiple refinement methods

**Sources**: pybader GitHub, PyPI documentation

## Key Strengths

### Python Integration:
- pip installable
- Scriptable workflows
- NumPy-based
- Jupyter compatible

### Performance:
- Multi-threaded execution
- Efficient memory usage
- Configurable precision
- Fast for routine analysis

### Flexibility:
- Multiple input formats
- Configurable methods
- Export options
- Easy automation

## Inputs & Outputs
- **Input formats**:
  - VASP CHGCAR
  - Cube files
  - Density grids
  
- **Output data types**:
  - Bader charges
  - Atomic volumes
  - Basin assignments
  - Charge density files

## Installation
```bash
pip install pybader
```

## Usage Examples
```python
from pybader import bader

# Run Bader analysis on CHGCAR
results = bader("CHGCAR")

# Access charges
charges = results.charges
volumes = results.volumes

# Command line usage
# pybader CHGCAR
```

## Performance Characteristics
- **Speed**: Competitive with compiled code
- **Memory**: Efficient grid handling
- **Parallelization**: Multi-threaded

## Limitations & Known Constraints
- **Grid-based**: Accuracy depends on grid density
- **Python overhead**: Slightly slower than Fortran
- **Large systems**: Memory for very large grids
- **Documentation**: Could be more extensive

## Comparison with Other Tools
- **vs Bader (Henkelman)**: pybader Python, Bader Fortran
- **vs Critic2**: pybader focused on Bader only
- **vs pymatgen Bader**: Different implementations
- **Unique strength**: Native Python, easy automation

## Application Areas
- Charge transfer analysis
- Oxidation state determination
- Bonding characterization
- High-throughput workflows
- Materials screening

## Best Practices
- Use fine charge density grids
- Include core charges for PAW
- Verify against known systems
- Check charge conservation

## Community and Support
- GitHub repository
- PyPI package
- MIT licensed
- Open development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/adam-kerrigan/pybader
2. PyPI: https://pypi.org/project/pybader/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- PyPI package: AVAILABLE
- Source code: OPEN (MIT)
- Method: Bader QTAIM analysis
