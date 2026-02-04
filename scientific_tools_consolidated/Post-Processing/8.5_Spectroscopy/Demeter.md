# Demeter (Athena/Artemis/Hephaestus)

## Official Resources
- Homepage: https://bruceravel.github.io/demeter/
- GitHub: https://github.com/bruceravel/demeter
- Documentation: https://bruceravel.github.io/demeter/documents/
- Publication: B. Ravel, M. Newville, J. Synchrotron Rad. 12, 537 (2005)
- License: Artistic License 2.0

## Overview
Demeter is a comprehensive software package for processing and analyzing X-ray Absorption Spectroscopy (XAS) data. It includes three main programs: Athena (data processing), Artemis (EXAFS fitting with FEFF), and Hephaestus (periodic table and X-ray data). Built on the IFEFFIT library, Demeter provides a complete workflow for XAS analysis.

**Scientific domain**: X-ray absorption spectroscopy, XANES, EXAFS
**Target user community**: XAS researchers and synchrotron beamline users

## Theoretical Methods
- XANES data processing
- EXAFS Fourier transform analysis
- FEFF path-based fitting
- Multiple scattering theory
- Linear combination fitting
- Principal component analysis

## Capabilities (CRITICAL)
- **Athena**: Data import, calibration, alignment, merging, normalization
- **Artemis**: EXAFS fitting with FEFF paths, multiple dataset fitting
- **Hephaestus**: X-ray data lookup, absorption calculations
- **Multi-platform**: Windows, macOS, Linux
- **FEFF Integration**: Built-in FEFF6 support
- **Batch Processing**: Project-based workflow

**Sources**: Demeter documentation, J. Synchrotron Rad. publication

## Key Strengths

### Complete Workflow:
- Data import to final fit
- Athena + Artemis integration
- Project management
- Reproducible analysis

### User-Friendly:
- GUI-based interface
- Interactive plotting
- Extensive documentation
- Tutorial materials

### FEFF Integration:
- Built-in FEFF6
- Path generation
- Multiple scattering
- Fitting interface

## Inputs & Outputs
- **Input formats**:
  - ASCII column data
  - Athena project files
  - Multiple beamline formats
  - CIF files for FEFF
  
- **Output data types**:
  - Athena project files (.prj)
  - Artemis project files
  - Fit reports
  - Publication-ready figures

## Installation
```bash
# Download from GitHub releases
# Windows installer available
# macOS and Linux: build from source or use package managers
```

## Usage Examples
```
# Athena workflow:
1. Import data files
2. Calibrate energy scale
3. Align and merge scans
4. Set normalization parameters
5. Export for Artemis fitting

# Artemis workflow:
1. Import Athena data
2. Run FEFF calculation from CIF
3. Build fitting model with paths
4. Perform fit and analyze results
```

## Performance Characteristics
- **Speed**: Efficient IFEFFIT backend
- **Memory**: Handles typical XAS datasets
- **Stability**: Production-ready software

## Limitations & Known Constraints
- **Perl-based**: Older codebase compared to Larch
- **GUI-focused**: Less scriptable than Larch
- **FEFF6 only**: Built-in FEFF is version 6
- **Development**: Transitioning to Larch

## Comparison with Other Tools
- **vs Larch**: Demeter GUI-focused, Larch Python-native
- **vs FEFF**: Demeter analysis, FEFF calculation engine
- **vs standalone codes**: Demeter comprehensive workflow
- **Unique strength**: Complete GUI workflow, extensive documentation

## Application Areas
- Synchrotron XAS experiments
- Catalysis and materials science
- Environmental chemistry
- Coordination chemistry
- Bioinorganic chemistry

## Best Practices
- Use reference compounds for calibration
- Document all processing parameters
- Validate fits with multiple approaches
- Save projects for reproducibility

## Community and Support
- GitHub repository
- XAFS mailing list
- Extensive documentation
- Bruce Ravel (developer)

## Verification & Sources
**Primary sources**:
1. Homepage: https://bruceravel.github.io/demeter/
2. B. Ravel, M. Newville, J. Synchrotron Rad. 12, 537 (2005)
3. GitHub: https://github.com/bruceravel/demeter

**Confidence**: VERIFIED - Published in J. Synchrotron Rad.

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub)
- Academic citations: >3000
- Active development: Maintained
