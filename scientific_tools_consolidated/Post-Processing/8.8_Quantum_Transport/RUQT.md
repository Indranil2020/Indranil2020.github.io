# RUQT

## Official Resources
- Source Repository: https://github.com/HoyLab-Rowan/RUQT
- Documentation: Included in repository
- License: Open source

## Overview
**RUQT** (Rowan University Quantum Transport) is a Fortran code for performing Landauer NEGF calculations using advanced electronic structure methods, particularly parametric 2-RDM (NEGF-RDM) and multi-configuration pair density functional theory (NEGF-MCPDFT) for quantum transport through molecular junctions.

**Scientific domain**: NEGF quantum transport, 2-RDM, MCPDFT for molecular junctions  
**Target user community**: Researchers studying quantum transport with advanced electronic structure methods beyond DFT

## Theoretical Methods
- Non-Equilibrium Green's Function (NEGF)
- Landauer formalism
- Parametric 2-RDM method (NEGF-RDM)
- Multi-configuration pair DFT (NEGF-MCPDFT)
- Molecular junction transport
- Beyond-DFT transport

## Capabilities (CRITICAL)
- NEGF quantum transport
- Landauer conductance
- 2-RDM transport (NEGF-RDM)
- MCPDFT transport (NEGF-MCPDFT)
- Molecular junction simulation
- Beyond-DFT accuracy

**Sources**: GitHub repository

## Key Strengths

### Beyond-DFT:
- 2-RDM method for transport
- MCPDFT for transport
- Correlation effects in transport
- Higher accuracy than DFT-NEGF

### Advanced Methods:
- NEGF-RDM framework
- NEGF-MCPDFT framework
- Multi-reference transport
- Strong correlation support

### Molecular Junctions:
- Molecular transport
- Junction I-V curves
- Transmission functions
- Contact effects

## Inputs & Outputs
- **Input formats**:
  - Molecular geometry
  - Electronic structure parameters
  - Junction configuration
  
- **Output data types**:
  - Transmission functions
  - I-V curves
  - Conductance
  - Current-voltage characteristics

## Interfaces & Ecosystem
- **Fortran**: Core computation
- **Python**: Wrapper scripts
- **Psi4/PySCF**: Electronic structure

## Performance Characteristics
- **Speed**: Moderate (2-RDM/MCPDFT)
- **Accuracy**: Beyond DFT
- **System size**: Small molecules
- **Memory**: Moderate

## Computational Cost
- **Transport**: Minutes to hours
- **Electronic structure**: Minutes (separate)
- **Typical**: Moderate

## Limitations & Known Constraints
- **Small molecules**: 2-RDM scaling limits size
- **Fortran compilation**: Required
- **Limited documentation**: Research code
- **Niche methods**: 2-RDM/MCPDFT community

## Comparison with Other Codes
- **vs Gollum**: RUQT has beyond-DFT methods, Gollum is DFT/TB-based
- **vs Transiesta**: RUQT uses 2-RDM/MCPDFT, Transiesta uses DFT
- **vs GreenCheetah**: RUQT is advanced electronic structure, GreenCheetah is recursive GF
- **Unique strength**: NEGF quantum transport with 2-RDM and MCPDFT beyond-DFT methods

## Application Areas

### Molecular Electronics:
- Molecular junction transport
- Single-molecule conductance
- I-V characteristics
- Contact effects

### Strongly Correlated Transport:
- Correlation effects in transport
- Multi-reference transport
- Transition metal complexes
- Open-shell systems

### Method Development:
- Beyond-DFT transport benchmarking
- 2-RDM transport validation
- MCPDFT transport testing
- Correlation-transport coupling

## Best Practices

### Electronic Structure:
- Use appropriate active space
- Check 2-RDM convergence
- Validate MCPDFT against experiment
- Compare with DFT-NEGF

### Transport:
- Use sufficient energy grid
- Check transmission convergence
- Validate I-V against experiment
- Consider contact effects

## Community and Support
- Open source on GitHub
- Developed at Rowan University (HoyLab)
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/HoyLab-Rowan/RUQT

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Specialized strength: NEGF quantum transport with 2-RDM and MCPDFT beyond-DFT methods
