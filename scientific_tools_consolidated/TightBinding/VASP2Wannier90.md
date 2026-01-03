# VASP2Wannier90

## Official Resources
- Homepage: Part of Wannier90 ecosystem
- Documentation: Wannier90 documentation, VASP wiki
- Source Repository: https://github.com/wannier-developers/wannier90 (interface components)
- License: Integrated with Wannier90 (GPL v2) and VASP (commercial)

## Overview
VASP2Wannier90 is the interface between the Vienna Ab initio Simulation Package (VASP) and Wannier90, enabling the construction of maximally-localized Wannier functions from VASP DFT calculations. This interface is one of the most widely used Wannier90 connections due to VASP's popularity, providing seamless workflow from VASP electronic structure calculations to Wannier tight-binding models.

**Scientific domain**: DFT-Wannier interface, VASP integration  
**Target user community**: VASP users, Wannier90 workflows, materials scientists

## Theoretical Methods
- DFT to Wannier transformation
- Overlap matrix calculation
- Projection construction
- MLWF generation
- Interface algorithms

## Capabilities (CRITICAL)
**Category**: DFT-Wannier interface
- VASP to Wannier90 interface
- Overlap matrix generation
- Projection calculation
- Wannier function construction
- Standard workflow component
- Production quality
- Widely used

**Sources**: Wannier90/VASP documentation

## Key Strengths

### VASP Integration:
- Native VASP support
- Standard workflow
- Production tested
- Widely used
- Well-documented

### Workflow:
- VASP DFT calculation
- Generate Wannier inputs
- Wannier90 processing
- Complete pipeline
- Automated

### Community Standard:
- Most common VASP-Wannier route
- Extensive validation
- Large user base
- Tutorial materials
- Best practices

## Workflow and Usage

### VASP INCAR Settings:
```
# Standard DFT calculation
LWANNIER90 = .TRUE.
LWRITE_MMN_AMN = .TRUE.

# For spin-orbit coupling
LSORBIT = .TRUE.
GGA_COMPAT = .FALSE.
```

### Wannier90 Preparation:
```bash
# Generate wannier90.win input
# Run Wannier90 in preprocessing mode
wannier90.x -pp silicon
```

### VASP Calculation:
```bash
# VASP generates .mmn, .amn files
mpirun -n 16 vasp_std
```

### Wannier90 Processing:
```bash
# Complete Wannierization
wannier90.x silicon
```

### Complete Workflow:
```bash
# 1. SCF calculation
# INCAR with standard DFT settings
vasp_std

# 2. NSCF on fine k-mesh
# INCAR with LWANNIER90 = .TRUE.
# Prepare wannier90.win
wannier90.x -pp silicon
vasp_std

# 3. Wannierize
wannier90.x silicon
```

## Advanced Features

### Spin-Orbit Coupling:
- LSORBIT support
- Spinor wavefunctions
- SOC Wannier functions
- Topological materials

### Hybrid Functionals:
- HSE, PBE0 support
- Accurate band structures
- Wannier from hybrid DFT
- High-quality MLWFs

### Projections:
- Atomic orbital projections
- Custom projections
- Disentanglement
- Optimal localization

## Performance Characteristics
- **Speed**: VASP-limited
- **Accuracy**: DFT quality
- **Integration**: Seamless
- **Purpose**: Production standard
- **Typical**: Standard DFT timescale

## Computational Cost
- VASP calculation dominant
- Interface overhead minimal
- Standard DFT workflow
- Production capable

## Limitations & Constraints
- **Requires VASP**: Commercial code
- **VASP version**: Compatibility requirements
- **Interface updates**: VASP version dependent
- **License**: VASP license needed

## Comparison with Other DFT-Wannier Interfaces
- **vs pw2wannier90 (QE)**: VASP2W90 for VASP, pw2w90 for QE
- **vs ABINIT**: VASP widely used commercially
- **Most common**: With VASP popularity
- **Standard**: Industry and academia

## Application Areas

### Materials Science:
- Electronic structure
- Tight-binding models
- Topological analysis
- Transport properties
- Universal application

### VASP Users:
- Standard post-processing
- Wannier workflows
- Property calculations
- Model construction

### Industrial:
- Commercial VASP users
- Production calculations
- Materials design
- Device modeling

## Best Practices

### VASP Settings:
- Appropriate k-mesh
- Converged DFT
- Proper INCAR flags
- VASP version compatibility

### Wannier90 Input:
- Proper projections
- Energy windows
- Disentanglement
- Validation

### Workflow:
- Test on small systems
- Check overlap matrices
- Validate spreads
- Compare band structures

## Community and Support
- VASP community
- Wannier90 community
- Combined documentation
- Large user base
- Forums and mailing lists

## Educational Resources
- VASP wiki
- Wannier90 tutorials
- Combined workflow guides
- Example calculations
- User forums

## Development
- VASP developers
- Wannier90 developers
- Joint interface maintenance
- Regular updates
- Version compatibility

## Research Impact
VASP2Wannier90 enables the vast VASP user community to construct Wannier functions, making Wannier-based analysis accessible to commercial and academic VASP users worldwide.

## Verification & Sources
**Primary sources**:
1. VASP wiki
2. Wannier90: https://wannier.org/
3. Interface documentation

**Secondary sources**:
1. User publications (thousands)
2. Tutorial materials

**Confidence**: VERIFIED - Standard VASP-Wannier interface

**Verification status**: ✅ VERIFIED
- **Category**: DFT-Wannier interface
- Status: Production standard
- **Note**: VASP2Wannier90 is the standard interface between VASP and Wannier90. Enables construction of maximally-localized Wannier functions from VASP DFT calculations. Most widely used Wannier90 interface due to VASP popularity. Production quality, well-documented, seamless workflow. Requires VASP license. Part of standard VASP-Wannier90 ecosystem.
