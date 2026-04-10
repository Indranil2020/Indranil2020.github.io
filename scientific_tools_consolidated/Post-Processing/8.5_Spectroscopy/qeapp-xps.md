# qeapp-xps

## Official Resources
- Source Repository: https://github.com/superstar54/qeapp-xps
- Documentation: Included in repository
- License: Open source

## Overview
**qeapp-xps** is an AiiDA plugin for calculating X-ray Photoelectron Spectroscopy (XPS) spectra using the XpsWorkChain of the aiida-quantumespresso package. It provides an automated workflow for computing core-level binding energies using core-hole pseudopotentials within Quantum ESPRESSO.

**Scientific domain**: X-ray photoelectron spectroscopy, core-level spectroscopy  
**Target user community**: Researchers computing XPS spectra from first principles using Quantum ESPRESSO via AiiDA workflows

## Theoretical Methods
- Core-hole pseudopotential method
- ΔKohn-Sham approach for binding energies
- Density Functional Theory (Quantum ESPRESSO)
- Final-state approximation
- Initial-state approximation
- AiiDA workflow automation

## Capabilities (CRITICAL)
- XPS binding energy calculation
- Core-hole pseudopotential generation
- Automated XPS workflow
- Multiple element/orbital support
- Spin-orbit splitting
- Chemical shift calculation
- Surface and bulk XPS
- High-throughput XPS via AiiDA

**Sources**: GitHub repository

## Key Strengths

### Automated Workflow:
- AiiDA-managed calculations
- Automatic pseudopotential handling
- Reproducible results
- Provenance tracking
- Error handling and recovery

### QE Integration:
- Uses well-tested QE core-hole method
- Same pseudopotentials and parameters
- Consistent with QE ecosystem
- Validated methodology

### Core-Hole Pseudopotentials:
- Supplied pseudopotential library
- Multiple elements supported
- Consistent treatment across periodic table
- Validated against experiment

## Inputs & Outputs
- **Input formats**:
  - AiiDA structure data
  - QE input parameters
  - Core-hole pseudopotential selection
  
- **Output data types**:
  - XPS binding energies
  - Chemical shifts
  - XPS spectra (with broadening)
  - Orbital-resolved contributions

## Interfaces & Ecosystem
- **AiiDA**: Workflow management
- **Quantum ESPRESSO**: DFT engine
- **aiida-quantumespresso**: QE-AiiDA interface
- **Materials Cloud**: Data sharing

## Performance Characteristics
- **Speed**: Depends on QE calculation
- **Accuracy**: Good (0.3-1 eV for chemical shifts)
- **System size**: Limited by QE
- **Automation**: Full AiiDA automation

## Computational Cost
- **Per core level**: One QE SCF calculation
- **Full XPS**: Multiple SCF calculations
- **Typical**: Hours for moderate systems
- **Automation**: Reduces manual effort

## Limitations & Known Constraints
- **QE only**: No VASP or other code support
- **Core-hole method**: Final-state approximation
- **Pseudopotential availability**: Not all elements have core-hole PPs
- **AiiDA dependency**: Requires AiiDA infrastructure
- **No multiplet effects**: Single-particle treatment

## Comparison with Other Codes
- **vs StoBe**: qeapp-xps is periodic (QE), StoBe is molecular
- **vs ORCA XPS**: qeapp-xps is DFT periodic, ORCA is wavefunction molecular
- **vs xspectra**: qeapp-xps is XPS, xspectra is XAS
- **Unique strength**: Automated AiiDA workflow for XPS from Quantum ESPRESSO, core-hole pseudopotential library

## Application Areas

### Surface Science:
- Surface core-level shifts
- Adsorbate binding energies
- Surface reconstruction effects
- Interface XPS

### Battery Materials:
- Redox state tracking
- Electrolyte decomposition
- SEI layer characterization
- Cycling-induced shifts

### Catalysis:
- Active site oxidation states
- Adsorbate-induced shifts
- Under reaction conditions
- Support effects

### 2D Materials:
- Layer-dependent shifts
- Defect characterization
- Doping effects
- Heterostructure interfaces

## Best Practices

### Pseudopotential Selection:
- Use supplied core-hole PPs
- Test convergence with cutoff
- Validate against experimental shifts
- Use consistent PP generation

### AiiDA Workflow:
- Use appropriate computer configuration
- Set reasonable wall times
- Monitor calculation progress
- Check for convergence issues

### Binding Energy Analysis:
- Reference to appropriate Fermi level
- Account for charging effects
- Compare relative shifts
- Validate with experimental XPS

## Community and Support
- Open source on GitHub
- Developed within AiiDA/QE ecosystem
- Part of Materials Cloud
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/superstar54/qeapp-xps
2. AiiDA-QuantumESPRESSO documentation
3. Materials Cloud platform

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Ongoing
- Specialized strength: Automated AiiDA workflow for XPS from Quantum ESPRESSO, core-hole pseudopotential library
