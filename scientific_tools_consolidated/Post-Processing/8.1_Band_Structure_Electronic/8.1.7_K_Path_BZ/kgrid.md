# kgrid

## Official Resources
- **Homepage**: https://github.com/WMD-group/kgrid
- **GitHub**: https://github.com/WMD-group/kgrid
- **PyPI**: https://pypi.org/project/kgrid/
- **License**: MIT License

## Overview
kgrid is a Python tool for calculating the required k-point density from input geometry for periodic quantum chemistry calculations. It uses a length cutoff approach to determine appropriate Monkhorst-Pack grids, ensuring consistent k-point density across different cell sizes and shapes.

**Scientific domain**: K-point sampling, DFT calculations, Brillouin zone integration
**Target user community**: DFT practitioners needing consistent k-point grids for convergence studies

## Theoretical Background
kgrid implements k-point selection based on:
- Real-space length cutoff: L_cutoff determines minimum sampling
- Monkhorst-Pack grid: N_i = ceil(L_cutoff / |a_i|)
- Reciprocal space density: Ensures consistent BZ sampling
- Accounts for cell shape and anisotropy

## Capabilities (CRITICAL)
- **K-point Generation**: Calculate optimal Monkhorst-Pack grids
- **Length Cutoff**: Specify real-space cutoff for consistent density
- **Multi-code Output**: VASP KSPACING, CASTEP MP_SPACING formats
- **Automatic Scaling**: Adjusts for cell dimensions

## Key Strengths

### Consistent Sampling:
- Length-based approach ensures comparable accuracy
- Automatic adjustment for cell shape
- Reproducible k-point selection

### Multi-Code Support:
- VASP KSPACING output
- CASTEP MP_SPACING format
- Easy integration into workflows

### Simple Interface:
- Command-line tool
- Python API
- Minimal dependencies

## Inputs & Outputs
- **Input formats**:
  - POSCAR (VASP structure)
  - Any ASE-readable structure
  - Length cutoff parameter
  
- **Output data types**:
  - KSPACING value (VASP)
  - MP_SPACING value (CASTEP)
  - k-point grid dimensions

## Installation
```bash
pip install kgrid
```

## Usage Examples
Command line:
```bash
# VASP format (default)
kgrid POSCAR 25  # 25 Å cutoff

# CASTEP format
kgrid --castep POSCAR 25

# Show grid dimensions
kgrid --verbose POSCAR 25
```

Python API:
```python
from kgrid import calc_kpt_tuple
from ase.io import read

atoms = read("POSCAR")
kpts = calc_kpt_tuple(atoms, cutoff_length=25)
print(f"K-point grid: {kpts}")
```

## Performance Characteristics
- **Speed**: Instant calculation
- **Accuracy**: Consistent density across systems
- **Simplicity**: Single parameter (cutoff length)

## Limitations & Known Constraints
- **Monkhorst-Pack only**: No support for other schemes
- **Gamma-centered**: Default behavior
- **Odd/even grids**: May need manual adjustment

## Comparison with Other Tools
- **vs SeeK-path**: kgrid for grids, SeeK-path for paths
- **vs manual selection**: kgrid ensures consistency
- **Unique strength**: Length-based consistent sampling

## Application Areas
- Convergence testing
- High-throughput calculations
- Consistent k-point selection
- Workflow automation

## Best Practices
- Use same cutoff for related calculations
- Verify convergence with increasing cutoff
- Consider symmetry for efficiency
- Document cutoff used for reproducibility

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/WMD-group/kgrid

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Developer: WMD Group (Walsh Materials Design)
