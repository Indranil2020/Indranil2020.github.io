# OCEAN (Obtaining Core Excitations using NBSE and Abinit)

## Official Resources
- Homepage: http://feff.phys.washington.edu/OCEAN/
- Documentation: http://feff.phys.washington.edu/OCEAN/documentation.html
- Source Repository: Available from developers
- License: Free for academic use

## Overview
OCEAN is a specialized code for calculating X-ray absorption spectra (XAS) and X-ray emission spectra (XES) using the Bethe-Salpeter Equation combined with DFT calculations. Developed at the University of Washington as part of the FEFF project, OCEAN focuses on core-level spectroscopy with emphasis on accurate treatment of core-hole effects, many-body interactions, and experimental comparison. It uses ABINIT or Quantum ESPRESSO for ground-state DFT and implements sophisticated BSE for core excitations.

**Scientific domain**: X-ray spectroscopy, core-level excitations, BSE  
**Target user community**: X-ray spectroscopists, synchrotron users, core-level researchers

## Theoretical Methods
- Bethe-Salpeter Equation (BSE)
- Core-level excitations
- DFT ground state (ABINIT/QE)
- Pseudopotentials and PAW
- Core-hole treatment
- Many-body effects
- X-ray absorption (XAS/XANES/NEXAFS)
- X-ray emission (XES)

## Capabilities (CRITICAL)
- X-ray absorption spectra (XAS)
- XANES (Near-edge structure)
- NEXAFS (Near-edge fine structure)
- X-ray emission spectra (XES)
- Core-level BSE
- Core-hole effects
- Element-specific spectra
- Polarization dependence
- Experimental comparison
- Production calculations

**Sources**: OCEAN website (http://feff.phys.washington.edu/OCEAN/)

## Key Strengths

### X-ray Spectroscopy Focus:
- Specialized for XAS/XES
- Core-level expertise
- Synchrotron applications
- Experimental validation
- Production quality

### Core-Hole Treatment:
- Explicit core-hole
- Many-body effects
- Accurate physics
- BSE formalism
- Sophisticated approach

### DFT Integration:
- ABINIT interface
- Quantum ESPRESSO interface
- Standard DFT input
- Flexible framework
- Wide compatibility

### FEFF Connection:
- Part of FEFF project
- Complementary tools
- Integrated workflow
- Expertise in X-ray
- Established lineage

### Experimental Comparison:
- Direct XAS comparison
- Synchrotron validation
- Quantitative predictions
- Edge calculations
- Practical applications

## Inputs & Outputs
- **Input formats**:
  - ABINIT/QE DFT output
  - OCEAN input files
  - Crystal structure
  - Core-hole specifications
  
- **Output data types**:
  - XAS spectra
  - XES spectra
  - Absorption edges
  - Polarization-dependent
  - Spectral functions

## Interfaces & Ecosystem
- **ABINIT**:
  - Primary DFT backend
  - PAW implementation
  - Tested workflow
  
- **Quantum ESPRESSO**:
  - Alternative DFT backend
  - Standard interface
  
- **FEFF**:
  - Related project
  - Complementary methods
  - Integrated ecosystem

## Workflow and Usage

### Typical Workflow:
1. Run DFT calculation (ABINIT/QE)
2. Prepare OCEAN input
3. Specify core-hole and edge
4. Run BSE calculation
5. Generate XAS/XES spectra
6. Compare with experiment

### Core-Level BSE:
- Core-hole excitation
- Electron-hole interaction
- Many-body screening
- Spectral calculation

## Advanced Features

### BSE Implementation:
- Core-level BSE
- Electron-hole kernel
- Many-body effects
- Accurate spectra
- Production quality

### Core-Hole Physics:
- Explicit core-hole
- Relaxation effects
- Screening
- Final-state interactions
- Comprehensive treatment

### Edge Calculations:
- K-edge, L-edge, M-edge
- Element-specific
- Polarization dependence
- Angular dependence
- Experimental geometry

### Spectroscopy Types:
- XAS (absorption)
- XANES (near-edge)
- NEXAFS (fine structure)
- XES (emission)
- Multiple techniques

## Performance Characteristics
- **Speed**: Moderate (BSE expense)
- **Accuracy**: Excellent for XAS
- **System size**: Moderate
- **Purpose**: Core spectroscopy
- **Typical**: Research and validation

## Computational Cost
- **BSE**: Expensive for core levels
- **DFT**: Standard cost
- **Production**: Feasible
- **Accuracy**: High
- **Value**: Experimental comparison

## Limitations & Known Constraints
- **Specialization**: Core spectroscopy focus
- **Computational cost**: BSE expense
- **DFT dependency**: Requires ABINIT/QE
- **Learning curve**: X-ray expertise helpful
- **Platform**: Linux systems

## Comparison with Other Codes
- **vs XSpectra (QE)**: OCEAN more sophisticated BSE
- **vs FEFF**: OCEAN BSE-based, FEFF multiple-scattering
- **vs Fiesta**: Both core-level, different emphases
- **Unique strength**: Specialized XAS/XES, core-hole BSE, synchrotron validation

## Application Areas

### Synchrotron Science:
- XAS measurements
- Beam line experiments
- Theoretical support
- Experimental interpretation
- Quantitative comparison

### Core-Level Spectroscopy:
- Element-specific
- Chemical sensitivity
- Oxidation states
- Local structure
- Electronic structure

### Materials Characterization:
- Electronic properties
- Chemical bonding
- Local environment
- Defects and impurities
- Surface science

## Best Practices

### DFT Preparation:
- Converged calculation
- PAW or appropriate pseudopotentials
- Core states treated properly
- Quality ground state

### Core-Hole Setup:
- Appropriate edge selection
- Core-hole specification
- Symmetry considerations
- Computational parameters

### Spectral Calculation:
- Convergence testing
- Broadening parameters
- Energy resolution
- Comparison with experiment

## Community and Support
- Free for academic use
- FEFF project support
- Documentation available
- University of Washington
- User community
- Synchrotron collaborations

## Educational Resources
- OCEAN documentation
- FEFF project materials
- XAS theory tutorials
- Example calculations
- Published papers

## Development
- University of Washington
- FEFF project team
- John Rehr group
- Active development
- Research-driven

## Research Applications
- X-ray absorption
- Core-level physics
- Synchrotron experiments
- Materials characterization
- Spectroscopy theory

## FEFF Project Integration
- Part of FEFF ecosystem
- Complementary to FEFF
- BSE approach vs multiple-scattering
- Integrated workflow
- Comprehensive X-ray tools

## Verification & Sources
**Primary sources**:
1. OCEAN website: http://feff.phys.washington.edu/OCEAN/
2. Documentation
3. University of Washington
4. FEFF project

**Secondary sources**:
1. X-ray spectroscopy literature
2. BSE for core levels
3. Synchrotron publications
4. Application studies

**Confidence**: VERIFIED - Established research code

**Verification status**: ✅ VERIFIED
- Website: ACCESSIBLE
- Documentation: Available
- Academic tool: Free for research
- FEFF project: CONFIRMED
- Development: University of Washington
- Specialized strength: X-ray absorption/emission spectroscopy, core-level BSE, core-hole treatment, synchrotron validation, ABINIT/QE integration, experimental comparison
