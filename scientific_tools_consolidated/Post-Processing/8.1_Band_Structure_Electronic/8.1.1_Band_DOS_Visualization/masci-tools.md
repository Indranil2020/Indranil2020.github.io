# masci-tools

## Official Resources
- **Homepage**: https://masci-tools.readthedocs.io/
- **GitHub**: https://github.com/JuDFTteam/masci-tools
- **Documentation**: https://masci-tools.readthedocs.io/
- **PyPI**: https://pypi.org/project/masci-tools/
- **License**: MIT License

## Overview
masci-tools is a Python package providing tools for pre- and post-processing of FLEUR and KKR (Korringa-Kohn-Rostoker) calculations. Developed by the JuDFT team at Forschungszentrum Jülich, it provides comprehensive support for all-electron full-potential methods.

**Scientific domain**: FLEUR/KKR post-processing, all-electron methods, magnetic materials
**Target user community**: FLEUR and KKR users, researchers studying magnetic and correlated systems

## Theoretical Background
masci-tools processes output from:
- FLEUR: Full-potential linearized augmented plane-wave (FLAPW)
- KKR: Korringa-Kohn-Rostoker Green's function method
- Magnetic properties: Spin-polarized and non-collinear magnetism
- Core-level spectroscopy

## Capabilities (CRITICAL)
- **FLEUR Support**: Full-potential LAPW pre/post-processing
- **KKR Support**: Green's function method analysis
- **Band Structure**: Electronic band plotting with character
- **DOS**: Total and projected density of states
- **Magnetic Properties**: Spin-resolved and non-collinear analysis
- **XML Handling**: Schema-based FLEUR input/output parsing
- **AiiDA Integration**: Workflow automation support

## Key Strengths

### FLEUR Integration:
- inp.xml parsing and generation
- out.xml result extraction
- Schema validation
- All FLEUR versions supported

### KKR Support:
- Green's function output analysis
- Impurity calculations
- Transport properties
- Magnetic exchange parameters

### AiiDA Workflows:
- aiida-fleur integration
- aiida-kkr integration
- Automated calculation chains

## Inputs & Outputs
- **Input formats**:
  - FLEUR inp.xml, out.xml
  - KKR output files
  - Structure files
  
- **Output data types**:
  - Band structures
  - DOS plots
  - Magnetic moments
  - Processed data arrays

## Installation
```bash
pip install masci-tools
```

## Usage Examples
```python
from masci_tools.io.parsers.fleur import outxml_parser
from masci_tools.vis.fleur import plot_fleur_bands

# Parse FLEUR output
data = outxml_parser("out.xml")

# Plot band structure
plot_fleur_bands("banddos.hdf", show=True)

# Extract magnetic moments
moments = data['magnetic_moments']
```

## Performance Characteristics
- **Speed**: Efficient XML parsing
- **Memory**: Handles large output files
- **Compatibility**: All FLEUR/KKR versions

## Limitations & Known Constraints
- **FLEUR/KKR specific**: Specialized for these codes
- **Schema dependency**: Requires matching XML schemas
- **Learning curve**: FLEUR/KKR knowledge needed

## Comparison with Other Tools
- **vs abipy**: masci-tools for FLEUR/KKR, abipy for ABINIT
- **vs pymatgen**: masci-tools specialized for all-electron
- **Unique strength**: FLEUR/KKR native support, AiiDA integration

## Application Areas
- All-electron calculations
- Magnetic materials
- Heavy elements (relativistic effects)
- Core-level spectroscopy
- Transport properties

## Verification & Sources
**Primary sources**:
1. Documentation: https://masci-tools.readthedocs.io/
2. GitHub: https://github.com/JuDFTteam/masci-tools

**Confidence**: VERIFIED - Official JuDFT team tool

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Developer: JuDFTteam (Forschungszentrum Jülich)
- Active development: Regular releases
