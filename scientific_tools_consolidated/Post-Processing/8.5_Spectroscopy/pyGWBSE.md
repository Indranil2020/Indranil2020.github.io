# pyGWBSE

## Official Resources
- Homepage: https://github.com/cmdlab/pyGWBSE
- GitHub: https://github.com/cmdlab/pyGWBSE
- Publication: Biswas et al., npj Comput. Mater. 9, 22 (2023)
- License: MIT License

## Overview
pyGWBSE is an open-source Python workflow package for performing high-throughput GW-BSE (Bethe-Salpeter equation) calculations. It automates the calculation of quasiparticle band structures, optical absorption spectra with excitonic effects, and effective masses using VASP as the DFT engine and interfaces with Wannier90 for band interpolation.

**Scientific domain**: GW-BSE calculations, optical spectra, excitons
**Target user community**: Researchers studying optical properties and excitons in materials

## Theoretical Methods
- G0W0 and GW0 quasiparticle calculations
- Bethe-Salpeter equation (BSE)
- Wannier interpolation for band structures
- Excitonic effects in optical spectra
- Effective mass calculation

## Capabilities (CRITICAL)
- **GW Calculations**: G0W0 and GW0 band structures
- **BSE Spectra**: Optical absorption with excitons
- **Wannier90 Integration**: Band interpolation
- **Effective Masses**: Electron and hole masses
- **High-Throughput**: Automated workflows
- **VASP Backend**: Uses VASP for DFT/GW/BSE
- **Convergence Tools**: Parameter testing

**Sources**: pyGWBSE GitHub, npj Comput. Mater. publication

## Key Strengths

### Automated Workflow:
- End-to-end automation
- Convergence testing
- Parameter optimization
- High-throughput ready

### VASP Integration:
- Direct VASP interface
- GW and BSE support
- Wannier90 coupling
- Standard workflow

### Analysis Tools:
- Band structure plotting
- Optical spectra
- Effective mass fitting
- Comparison utilities

## Inputs & Outputs
- **Input formats**:
  - VASP input files
  - Materials Project IDs
  - pymatgen structures
  
- **Output data types**:
  - Quasiparticle band structures
  - Optical absorption spectra
  - Effective masses
  - JSON/HDF5 data

## Installation
```bash
pip install pyGWBSE
# Requires VASP and Wannier90
```

## Usage Examples
```python
from pyGWBSE import GWBSE

# Initialize workflow
calc = GWBSE(
    structure='POSCAR',
    functional='PBE',
    gw_type='G0W0'
)

# Run GW calculation
calc.run_gw()

# Run BSE for optical spectra
calc.run_bse()

# Get absorption spectrum
spectrum = calc.get_absorption()
```

## Performance Characteristics
- **Speed**: Depends on VASP GW/BSE
- **Memory**: GW/BSE memory intensive
- **Scalability**: High-throughput capable

## Limitations & Known Constraints
- **VASP required**: Commercial DFT code needed
- **Computational cost**: GW/BSE expensive
- **Convergence**: Careful parameter testing needed
- **Memory**: Large memory for GW

## Comparison with Other Tools
- **vs Yambo**: pyGWBSE VASP-based, Yambo QE-based
- **vs BerkeleyGW**: Different backends
- **vs manual workflow**: pyGWBSE automated
- **Unique strength**: High-throughput VASP GW-BSE

## Application Areas
- Semiconductor optical properties
- 2D material excitons
- Solar cell materials
- Optoelectronic devices
- Band gap engineering

## Best Practices
- Converge k-points carefully
- Test NBANDS for GW
- Validate against experimental spectra
- Document convergence parameters

## Community and Support
- GitHub repository
- npj Comput. Mater. publication
- MIT licensed
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/cmdlab/pyGWBSE
2. Biswas et al., npj Comput. Mater. 9, 22 (2023)

**Confidence**: VERIFIED - Published in npj Computational Materials

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Publication: npj Comput. Mater. (2023)
- Source code: OPEN (MIT)
- Method: GW-BSE workflow
