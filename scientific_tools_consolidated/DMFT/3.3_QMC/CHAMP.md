# CHAMP (Cornell-Holland Ab-initio Materials Package)

## Official Resources
- Homepage: https://github.com/CHAMPlib/CHAMP
- Documentation: GitHub repository and wiki
- Source Repository: https://github.com/CHAMPlib/CHAMP
- License: GNU General Public License v3.0

## Overview
CHAMP is a quantum Monte Carlo package originally developed at Cornell University and currently maintained as an open-source community project. The code implements Variational Monte Carlo (VMC) and Diffusion Monte Carlo (DMC) methods for electronic structure calculations of molecules and solids. CHAMP emphasizes flexibility, ease of modification for research, and educational value, making it suitable for both production calculations and QMC method development.

**Scientific domain**: Quantum Monte Carlo, electronic structure, ab-initio calculations  
**Target user community**: QMC researchers, method developers, educational users

## Theoretical Methods
- Variational Monte Carlo (VMC)
- Diffusion Monte Carlo (DMC)
- Slater-Jastrow wavefunctions
- Multi-determinant trial functions
- Pseudopotentials
- All-electron calculations
- Fixed-node approximation
- Wavefunction optimization

## Capabilities (CRITICAL)
**Category**: Open-source QMC code
- VMC and DMC methods
- Molecules and solids
- Periodic boundary conditions
- Slater-Jastrow wavefunctions
- Multi-determinant expansions
- Wavefunction optimization
- Energy calculations
- Forces
- Pseudopotentials
- Finite/periodic systems
- Research-friendly code
- Educational applications

**Sources**: GitHub repository, documentation

## Key Strengths

### Research-Friendly:
- Clear code structure
- Easy modification
- Method development
- Educational value
- Open-source

### Community Code:
- GitHub-based
- Community contributions
- Active development
- Issue tracking
- Collaborative

### Flexibility:
- Various trial functions
- Multiple implementations
- Research extensions
- Custom features
- Development platform

## Inputs & Outputs
- **Input formats**:
  - CHAMP input files
  - DFT trial wavefunctions (various)
  - Pseudopotentials
  - Structure files
  
- **Output data types**:
  - Total energies
  - Forces
  - Observables
  - Statistical data
  - Wavefunction parameters

## Interfaces & Ecosystem

### DFT Codes:
- GAMESS
- Gaussian
- PySCF
- Various converters

### Community:
- CHAMPlib organization
- GitHub collaboration
- User contributions
- Development community

## Workflow and Usage

### Installation:
```bash
# Clone repository
git clone https://github.com/CHAMPlib/CHAMP.git
cd CHAMP
# Follow build instructions
make
```

### Basic VMC:
```bash
# Prepare input files
# Run CHAMP
champ < input.inp > output.out
```

### DMC Calculation:
```bash
# Setup DMC input
# Run calculation
champ_dmc < dmc.inp > dmc.out
```

## Advanced Features

### Wavefunction Types:
- Slater-Jastrow
- Multi-determinant
- Various Jastrow forms
- Optimization tools

### Forces:
- DMC forces
- Geometry optimization
- Structural properties

### Method Development:
- Research modifications
- Algorithm testing
- Custom features
- Development platform

## Performance Characteristics
- **Speed**: Moderate (research code)
- **Accuracy**: QMC quality
- **Purpose**: Research and education
- **Scalability**: MPI support

## Computational Cost
- Standard QMC cost
- Research applications
- Production capable
- Educational use

## Limitations & Known Constraints
- **Performance**: Not most optimized
- **Documentation**: GitHub-based
- **Community**: Smaller than QMCPACK/CASINO
- **HPC optimization**: Less than major codes
- **Best for**: Research, education, development

## Comparison with Other QMC Codes
- **vs QMCPACK**: CHAMP research-friendly, QMCPACK production
- **vs CASINO**: CHAMP open/modifiable, CASINO feature-rich
- **Unique strength**: Open development, educational value, research flexibility, community code

## Application Areas

### Research:
- Method development
- Algorithm testing
- QMC research
- Custom features

### Education:
- Learning QMC
- Teaching tool
- Code understanding
- Student projects

### Production:
- Small to medium calculations
- Standard QMC
- Validation studies

## Best Practices

### Usage:
- Understand code structure
- Start with examples
- Community support
- GitHub issues

### Development:
- Contribute back
- Code documentation
- Testing
- Collaboration

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- CHAMPlib organization
- Issue tracking
- Community-driven
- Educational focus

## Educational Resources
- GitHub wiki
- Example inputs
- Source code (educational)
- User contributions
- QMC literature

## Development
- Community project
- Cornell origin
- Open development
- Active contributions
- Research focus

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/CHAMPlib/CHAMP
2. Repository documentation

**Secondary sources**:
1. QMC literature
2. User publications

**Confidence**: VERIFIED - Community QMC code

**Verification status**: ✅ VERIFIED
- GitHub: ACCESSIBLE
- License: GPL v3 (open-source)
- **Category**: Open-source QMC code
- Status: Community maintained
- Specialized strength: Research-friendly quantum Monte Carlo, VMC/DMC methods, educational value, open development, community code, method development platform, flexible for modifications, Cornell origin, GitHub-based collaboration
