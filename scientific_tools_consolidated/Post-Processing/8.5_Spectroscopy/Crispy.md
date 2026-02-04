# Crispy

## Official Resources
- Homepage: https://www.esrf.fr/computing/scientific/crispy/
- GitHub: https://github.com/mretegan/crispy
- Documentation: https://www.esrf.fr/computing/scientific/crispy/
- Publication: M. Retegan, Crispy: v0.7.4 (2024)
- License: MIT License

## Overview
Crispy is a modern graphical user interface for calculating core-level spectra using the semi-empirical multiplet approaches implemented in Quanty. It enables simulation of XAS (X-ray Absorption Spectroscopy), XES (X-ray Emission Spectroscopy), RIXS (Resonant Inelastic X-ray Scattering), and XPS spectra for transition metal and rare earth compounds.

**Scientific domain**: Core-level spectroscopy simulation (XAS, XES, RIXS)
**Target user community**: X-ray spectroscopists at synchrotrons and labs

## Theoretical Methods
- Atomic multiplet theory (via Quanty)
- Crystal field effects
- Charge transfer multiplets
- Spin-orbit coupling
- Hybridization effects
- RIXS cross-sections

## Capabilities (CRITICAL)
- **XAS Simulation**: L-edge, M-edge, K-edge
- **XES Simulation**: X-ray emission spectra
- **RIXS Maps**: 2D RIXS plane calculations
- **XPS Simulation**: Core-level photoemission
- **GUI Interface**: User-friendly Qt interface
- **Quanty Backend**: Full multiplet calculations
- **Cross-Platform**: Windows, macOS, Linux

**Sources**: Crispy documentation, ESRF development

## Key Strengths

### Modern GUI:
- Qt-based interface
- Interactive parameters
- Real-time preview
- Export capabilities

### Quanty Integration:
- Full multiplet theory
- Charge transfer
- RIXS support
- Validated calculations

### Cross-Platform:
- Windows, macOS, Linux
- Standalone installers
- Python package
- Active development

## Inputs & Outputs
- **Input formats**:
  - GUI parameter entry
  - Quanty input files
  - Python scripting
  
- **Output data types**:
  - Simulated spectra
  - RIXS maps
  - ASCII data files
  - Publication figures

## Installation
```bash
pip install crispy
# Or download standalone installer from GitHub
```

## Usage Examples
```python
# GUI usage:
# 1. Launch Crispy
# 2. Select element (e.g., Ni)
# 3. Choose edge (e.g., L2,3)
# 4. Set symmetry and crystal field
# 5. Run calculation
# 6. Compare with experiment

# Python API also available
from crispy import Calculation
calc = Calculation(element='Ni', edge='L2,3')
calc.run()
```

## Performance Characteristics
- **Speed**: Fast Quanty backend
- **Memory**: Depends on calculation size
- **Usability**: Intuitive GUI

## Limitations & Known Constraints
- **Quanty required**: Backend dependency
- **Semi-empirical**: Parameter fitting needed
- **Cluster size**: Limited to small clusters
- **Learning curve**: Multiplet concepts required

## Comparison with Other Tools
- **vs CTM4XAS**: Crispy modern GUI, more features
- **vs EDRIXS**: Different theoretical approaches
- **vs raw Quanty**: Crispy user-friendly wrapper
- **Unique strength**: Modern GUI for Quanty calculations

## Application Areas
- Synchrotron XAS beamlines
- RIXS experiments
- Transition metal chemistry
- Rare earth compounds
- Magnetic materials

## Best Practices
- Start with literature parameters
- Validate XAS before RIXS
- Consider charge transfer effects
- Compare different symmetries

## Community and Support
- GitHub repository
- ESRF development
- MIT licensed
- Marius Retegan (developer)

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mretegan/crispy
2. Documentation: https://www.esrf.fr/computing/scientific/crispy/

**Confidence**: VERIFIED - ESRF developed

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (MIT)
- Developer: ESRF
- Active development: Maintained
