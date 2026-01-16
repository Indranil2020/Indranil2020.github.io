# eT (electronic Transport)

## Official Resources
- Homepage: https://etprogram.org/
- Documentation: https://etprogram.org/user_manual.html
- Source Repository: https://gitlab.com/eT-program/eT
- License: GNU General Public License v3.0

## Overview
eT is a quantum chemistry program specialized in calculating molecular response properties and time-dependent phenomena using coupled cluster theory. Developed with focus on response properties, excited states, and spectroscopic calculations, eT implements efficient coupled cluster methods with particular emphasis on equations-of-motion and linear response approaches for accurate molecular properties.

**Scientific domain**: Coupled cluster, response properties, excited states, spectroscopy  
**Target user community**: Spectroscopy researchers, response property calculations, CC method users

## Theoretical Methods
- Coupled cluster (CC2, CCSD, CC3)
- Equations-of-motion coupled cluster (EOM-CC)
- Linear response coupled cluster
- Time-dependent coupled cluster
- Hartree-Fock
- Response properties
- Excited states
- Transition properties
- Spectroscopic properties

## Capabilities (CRITICAL)
- Ground-state coupled cluster
- Excited states (EOM-CC)
- Response properties
- Molecular properties
- Spectroscopic calculations
- Transition moments
- Oscillator strengths
- UV/Vis spectra
- Linear response
- Time-dependent properties
- High-accuracy calculations
- Research-quality results

**Sources**: Official website (https://etprogram.org/)

## Key Strengths

### Response Properties:
- Specialized focus
- CC response theory
- Accurate properties
- Transition moments
- Spectroscopic data

### Coupled Cluster:
- CC2, CCSD, CC3
- EOM-CC for excitations
- High accuracy
- Benchmark quality
- Production methods

### Spectroscopy:
- UV/Vis calculations
- Excitation energies
- Oscillator strengths
- Transition properties
- Experimental comparison

### Open Source:
- GPL v3 licensed
- GitLab repository
- Free to use
- Community development
- Transparent code

### Research Tool:
- Method development
- Property calculations
- Benchmark studies
- Accurate results
- Academic focus

## Inputs & Outputs
- **Input formats**:
  - Text-based input
  - Molecular coordinates
  - Basis set specifications
  - Method settings
  
- **Output data types**:
  - Energies
  - Excitation energies
  - Properties
  - Transition moments
  - Spectroscopic data

## Interfaces & Ecosystem
- **Development**:
  - GitLab repository
  - Community contributions
  - Research tool
  
- **Analysis**:
  - Property extraction
  - Spectroscopy analysis
  - Standard tools
  - Custom scripts

## Workflow and Usage

### Typical Workflow:
1. Prepare molecular structure
2. Select CC method
3. Choose basis set
4. Configure response properties
5. Run eT calculation
6. Analyze results

### Running eT:
```bash
et input.inp
# Runs coupled cluster calculation
```

## Advanced Features

### EOM-CC:
- Equations-of-motion
- Excited states
- Accurate excitations
- Multiple states
- Transition properties

### Response Theory:
- Linear response
- Time-dependent properties
- Molecular properties
- Spectroscopic calculations
- High accuracy

### CC Methods:
- CC2 (approximate)
- CCSD (standard)
- CC3 (high accuracy)
- Gradients (selected)
- Production quality

## Performance Characteristics
- **Speed**: Moderate (high-level CC)
- **Accuracy**: Excellent for properties
- **System size**: Small to medium molecules
- **Scalability**: Typical CC scaling
- **Typical**: Research calculations

## Computational Cost
- **CC2**: Moderate cost
- **CCSD**: Expensive
- **CC3**: Very expensive
- **Response**: Additional cost
- **Production**: Research-level feasible

## Limitations & Known Constraints
- **System size**: Limited by CC scaling
- **Community**: Smaller, specialized
- **Documentation**: Research-level
- **Features**: Focused on properties
- **Platform**: Linux primarily
- **Development**: Academic/research

## Comparison with Other Codes
- **vs CFOUR**: eT specialized for response properties
- **vs Dalton**: Similar focus, different implementation
- **vs PSI4**: eT more specialized
- **Unique strength**: CC response properties, EOM-CC, spectroscopy focus

## Application Areas

### Spectroscopy:
- UV/Vis calculations
- Excitation energies
- Transition moments
- Oscillator strengths
- Experimental comparison

### Response Properties:
- Molecular properties
- Linear response
- Time-dependent
- High accuracy
- Benchmark calculations

### Excited States:
- EOM-CC calculations
- Multiple states
- Accurate excitations
- State properties
- Transition densities

### Method Development:
- CC response research
- Algorithm testing
- Benchmark studies
- Property methods

## Best Practices

### Method Selection:
- CC2 for screening
- CCSD for production
- CC3 for benchmarks
- Appropriate for system
- Balance accuracy/cost

### Basis Sets:
- Appropriate for properties
- Augmented when needed
- Convergence testing
- Balance size/accuracy

### Response Calculations:
- Select properties carefully
- Converge SCF well
- Check stability
- Validate results

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Academic development
- Research community
- Limited production support

## Educational Resources
- GitHub documentation
- Source code
- Published papers
- CC theory literature
- Response property references

## Development
- GitHub-based development
- Academic contributors
- Research focus
- Ongoing improvements
- Community input

## Research Applications
- Response property calculations
- Spectroscopy predictions
- Benchmark studies
- Method development
- High-accuracy properties

## Technical Innovation

### CC Response:
- Efficient implementations
- Response theory
- Property calculations
- Specialized algorithms
- Research quality

### EOM-CC:
- Excited state methods
- Accurate excitations
- Transition properties
- Multiple states
- Production implementations

## Verification & Sources
**Primary sources**:
1. Official website: https://etprogram.org/
2. Documentation: https://etprogram.org/user_manual.html
3. GitLab repository: https://gitlab.com/eT-program/eT
3. Published papers on eT methodology

**Secondary sources**:
1. Coupled cluster literature
2. Response property methods
3. EOM-CC references
4. Spectroscopy calculations

**Confidence**: LOW_CONF - Research code, specialized focus, smaller community

**Verification status**: ✅ VERIFIED
- Official website: ACCESSIBLE (https://etprogram.org/)
- Documentation: ACCESSIBLE (user manual)
- Source code: OPEN (GitLab, GPL v3)
- Community support: GitLab, academic
- Development: Active on GitLab
- Specialized strength: Coupled cluster response properties, EOM-CC excited states, spectroscopic calculations, molecular properties, linear response theory, benchmark-quality accuracy
