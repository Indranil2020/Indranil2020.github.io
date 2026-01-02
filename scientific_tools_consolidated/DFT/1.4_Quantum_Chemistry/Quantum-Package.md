# Quantum Package

## Official Resources
- Homepage: https://quantumpackage.github.io/qp2/
- Documentation: https://quantumpackage.github.io/qp2/
- Source Repository: https://github.com/QuantumPackage/qp2
- License: GNU Affero General Public License v3.0

## Overview
Quantum Package is a programming environment for quantum chemistry developed in France, focusing on wavefunction methods and quantum Monte Carlo. Developed primarily at the Laboratoire de Chimie et Physique Quantiques in Toulouse, Quantum Package provides a modular framework for implementing and developing wavefunction-based methods with emphasis on selected configuration interaction and stochastic approaches.

**Scientific domain**: Wavefunction methods, selected CI, quantum Monte Carlo, method development  
**Target user community**: Method developers, wavefunction theory researchers, QMC specialists

## Theoretical Methods
- Selected configuration interaction (SCI)
- CIPSI (Configuration Interaction using a Perturbative Selection made Iteratively)
- Full configuration interaction (FCI)
- Quantum Monte Carlo (VMC, DMC)
- Hartree-Fock
- Density Functional Theory
- Coupled cluster
- Multi-reference methods
- Stochastic approaches
- Determinant-based methods

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Selected CI calculations
- CIPSI method
- Full CI for small systems
- Quantum Monte Carlo
- Excited states
- Multi-reference wavefunctions
- Benchmark-quality accuracy
- Stochastic methods
- Method development platform
- Modular architecture
- Python interface
- Parallel execution
- Research and development tool

**Sources**: GitHub repository (https://github.com/QuantumPackage/qp2)

## Key Strengths

### Selected CI:
- CIPSI algorithm
- Efficient wavefunction compression
- Controlled accuracy
- Large active spaces
- Systematic convergence

### Modularity:
- Plugin architecture
- Easy method development
- Extensible framework
- Research platform
- Custom modules

### QMC:
- Variational Monte Carlo
- Diffusion Monte Carlo
- Stochastic approaches
- High accuracy
- Benchmark quality

### Open Source:
- AGPL v3 licensed
- GitHub repository
- Free to use
- Community development
- Transparent code

### Method Development:
- Research tool
- Algorithm testing
- New method implementation
- Flexible framework
- Educational value

## Inputs & Outputs
- **Input formats**:
  - EZFIO database format
  - XYZ coordinates
  - Basis set specifications
  - Python interface
  
- **Output data types**:
  - Energies and wavefunctions
  - CI coefficients
  - Molecular orbitals
  - QMC data
  - Analysis files

## Interfaces & Ecosystem
- **Python Interface**:
  - Python API
  - Scripting capabilities
  - Workflow automation
  - Custom analysis
  
- **EZFIO**:
  - Hierarchical database
  - Data storage
  - Input/output management
  - Efficient access
  
- **Development**:
  - Plugin system
  - Module development
  - GitHub collaboration
  - Community contributions

## Workflow and Usage

### Typical Workflow:
1. Create EZFIO database
2. Run SCF calculation
3. Select CI method (CIPSI)
4. Run wavefunction calculation
5. Analyze results
6. Optional: QMC for refinement

### Command-Line Usage:
```bash
# Create database
qp_create_ezfio molecule.xyz -b cc-pvdz

# Run Hartree-Fock
qp_run scf molecule.ezfio

# Run CIPSI
qp_run fci molecule.ezfio
```

### Python Interface:
```python
from qp import *
# Python scripting for workflows
```

## Advanced Features

### CIPSI:
- Iterative selection
- Perturbative correction
- Systematic improvement
- Efficient for large systems
- Benchmark quality

### Stochastic Methods:
- Monte Carlo sampling
- Stochastic CI
- Variance reduction
- Controlled accuracy
- Scalable

### Multi-Reference:
- Large active spaces
- Multiple configurations
- Strongly correlated systems
- Accurate description
- Flexible selection

### Plugin System:
- Custom modules
- Method development
- Research extensions
- Community plugins
- Easy integration

## Performance Characteristics
- **Speed**: Varies by method
- **Accuracy**: Benchmark quality for selected methods
- **System size**: Small to medium molecules
- **Scalability**: Good parallelization
- **Typical**: Research calculations

## Computational Cost
- **CIPSI**: Expensive but efficient
- **FCI**: Very expensive (small systems)
- **QMC**: Stochastic, time-dependent
- **Development focus**: Not production speed
- **Research**: Accuracy-focused

## Limitations & Known Constraints
- **Speed**: Not optimized for production
- **Documentation**: Research-level
- **Learning curve**: Steep
- **Community**: Smaller, specialized
- **Focus**: Method development vs production
- **Platform**: Linux primarily
- **Maturity**: Research tool

## Comparison with Other Codes
- **vs Production codes**: QP development-focused
- **vs MOLPRO/MOLCAS**: QP more modular, research-oriented
- **vs GAMESS/NWChem**: QP specialized for CI/QMC development
- **Unique strength**: Selected CI (CIPSI), modularity, QMC, method development platform

## Application Areas

### Method Development:
- Algorithm research
- New methods
- Testing approaches
- Benchmarking
- Code prototyping

### Benchmark Calculations:
- Reference energies
- Accuracy standards
- Method validation
- Small system exactness
- Quality assessment

### Strongly Correlated:
- Multi-reference systems
- Bond breaking
- Transition metals
- Selected CI advantages
- Accurate treatment

### QMC Studies:
- Quantum Monte Carlo
- Stochastic methods
- High-accuracy energies
- Research applications

## Best Practices

### Method Selection:
- CIPSI for general use
- FCI for small systems
- QMC for refinement
- Appropriate for system
- Systematic approach

### Convergence:
- CIPSI threshold
- Selection criteria
- PT2 energy
- Systematic improvement
- Check convergence

### Development:
- Use plugin system
- Follow architecture
- Contribute to GitHub
- Test thoroughly
- Document code

## Community and Support
- Open-source (AGPL v3)
- GitHub repository
- French research community
- User mailing list
- Active development
- Research collaborations

## Educational Resources
- GitHub documentation
- Wiki pages
- Example calculations
- Published papers
- Tutorials (growing)

## Development
- Laboratoire de Chimie et Physique Quantiques (Toulouse)
- Anthony Scemama
- Michel Caffarel
- French research groups
- International collaborations
- Active GitHub development

## Research Focus

### CIPSI Method:
- Selected configuration interaction
- Perturbative selection
- Iterative improvement
- Benchmark quality
- Research applications

### Stochastic Approaches:
- Monte Carlo methods
- Variance reduction
- Scalable algorithms
- High accuracy
- Active research

### Modularity:
- Framework for development
- Easy prototyping
- Method testing
- Community contributions
- Extensible design

## French Development
- Strong French tradition
- Toulouse expertise
- European collaboration
- Academic excellence
- Research-driven

## Technical Innovation

### Selected CI:
- Efficient wavefunction representation
- Systematic selection
- Perturbative corrections
- Scalable approach
- Benchmark quality

### Modular Architecture:
- Plugin-based
- EZFIO database
- Clean separation
- Easy extension
- Research platform

## Verification & Sources
**Primary sources**:
1. Homepage: https://quantumpackage.github.io/qp2/
2. GitHub: https://github.com/QuantumPackage/qp2
3. A. Scemama et al., J. Chem. Theory Comput. papers on Quantum Package
4. CIPSI method publications

**Secondary sources**:
1. GitHub documentation
2. Published studies using Quantum Package
3. Selected CI literature
4. French quantum chemistry community

**Confidence**: LOW_CONF - Research tool, specialized methods, smaller community

**Verification status**: ✅ VERIFIED
- GitHub: ACCESSIBLE
- Documentation: Available online
- Source code: OPEN (GitHub, AGPL v3)
- Community support: Mailing list, GitHub
- Academic citations: Growing
- Active development: Regular GitHub activity
- Specialized strength: Selected CI (CIPSI), quantum Monte Carlo, modular architecture, method development platform, wavefunction methods, research tool, stochastic approaches
