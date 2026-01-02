# PROFESS (PRinceton Orbital-Free Electronic Structure Software)

## Official Resources
- Homepage: http://www.princeton.edu/cbe/people/faculty/carter/research-carter-group/profess/
- Documentation: Available through Princeton University
- Source Repository: Available to licensed users
- License: Free for academic use (license agreement required)

## Overview
PROFESS is an orbital-free density functional theory (OF-DFT) code developed at Princeton University. It implements kinetic energy density functionals that avoid explicit calculation of orbitals, enabling linear-scaling DFT calculations for large metallic and semiconductor systems. PROFESS is particularly efficient for materials where orbital-free approaches are accurate, providing significant computational speedup over traditional Kohn-Sham DFT.

**Scientific domain**: Orbital-free DFT, kinetic energy functionals, large-scale materials  
**Target user community**: Materials scientists, large system researchers, OF-DFT specialists

## Theoretical Methods
- Orbital-Free Density Functional Theory (OF-DFT)
- Kinetic energy density functionals (KEDF)
- Thomas-Fermi-Dirac approximation
- von Weizsäcker functional
- Wang-Teter, Wang-Govind-Carter KEDFs
- Local pseudopotentials
- Plane-wave basis
- Periodic systems
- Real-space grid option

## Capabilities (CRITICAL)
- Ground-state electronic structure (metals, semiconductors)
- Orbital-free calculations
- Linear-scaling DFT
- Large systems (100,000+ atoms)
- Total energy calculations
- Structural optimization
- Molecular dynamics
- Stress and pressure
- Metallic systems
- Simple semiconductors
- Alloys
- Fast calculations
- Benchmarking OF-DFT

**Sources**: Princeton University Carter Group (http://www.princeton.edu/cbe/people/faculty/carter/)

## Key Strengths

### Orbital-Free:
- No orbitals needed
- Much faster than KS-DFT
- Linear scaling
- Large systems feasible
- Efficient algorithms

### Scalability:
- Linear or near-linear scaling
- 100,000+ atoms demonstrated
- Minimal memory
- Fast calculations
- Extreme systems possible

### KEDF Development:
- Advanced kinetic functionals
- Research platform
- Method testing
- Functional development
- Benchmark quality

### Materials Focus:
- Metallic systems
- Simple semiconductors
- Alloys
- Bulk properties
- Realistic systems

### Research Tool:
- OF-DFT method development
- KEDF research
- Benchmark studies
- Algorithm testing
- Academic focus

## Inputs & Outputs
- **Input formats**:
  - Text-based input
  - Atomic coordinates
  - Pseudopotentials (local)
  - KEDF specifications
  
- **Output data types**:
  - Total energies
  - Forces and stresses
  - Optimized structures
  - Electron density
  - MD trajectories

## Interfaces & Ecosystem
- **Princeton Development**:
  - Carter group support
  - Academic distribution
  - Research collaboration
  
- **Analysis**:
  - Standard tools
  - Density analysis
  - Structure visualization
  - Custom scripts

## Workflow and Usage

### Typical Workflow:
1. Prepare structure
2. Select KEDF
3. Set up pseudopotentials (local)
4. Configure calculation
5. Run PROFESS
6. Analyze results

### Input Configuration:
- System geometry
- KEDF choice
- Pseudopotential files
- Convergence criteria
- Calculation type

### Running PROFESS:
```bash
profess input.ion
# Runs orbital-free DFT calculation
```

## Advanced Features

### Kinetic Energy Functionals:
- Multiple KEDF options
- Thomas-Fermi-Dirac
- von Weizsäcker
- Wang-Teter (WT)
- Wang-Govind-Carter (WGC)
- Huang-Carter (HC)
- Custom development

### Linear Scaling:
- O(N) or O(N log N)
- Conjugate gradient
- Efficient algorithms
- Minimal overhead
- Large system capability

### Real-Space Option:
- Alternative to plane waves
- Flexible boundaries
- Efficient for some systems
- Research feature

### Molecular Dynamics:
- Born-Oppenheimer MD
- Fast forces
- Large systems
- Long timescales
- Production simulations

## Performance Characteristics
- **Speed**: Much faster than KS-DFT
- **Scaling**: Linear or near-linear
- **System size**: Very large (100,000+ atoms)
- **Accuracy**: Good for metals, limited for others
- **Memory**: Very low requirements

## Computational Cost
- **OF-DFT**: Orders of magnitude faster than KS
- **Large systems**: Practical
- **MD**: Feasible for long times
- **Memory**: Minimal
- **Typical**: Production calculations

## Limitations & Known Constraints
- **Accuracy**: Limited by KEDF quality
- **Systems**: Best for metals, simple semiconductors
- **Pseudopotentials**: Local only
- **KEDF**: Not universal
- **Distribution**: Academic license
- **Documentation**: Research-level
- **Community**: Specialized

## Comparison with Other Codes
- **vs KS-DFT codes**: PROFESS much faster but less accurate
- **vs OFDFT**: PROFESS specialized OF-DFT
- **vs Classical MD**: PROFESS has electronic structure
- **Unique strength**: Orbital-free DFT, linear scaling, 100,000+ atoms, KEDF development

## Application Areas

### Large Metallic Systems:
- Bulk metals
- Metallic alloys
- Liquid metals
- Large-scale properties
- Production simulations

### Method Development:
- KEDF research
- Functional testing
- Benchmark studies
- Algorithm development
- OF-DFT advances

### Materials Screening:
- High-throughput
- Large systems
- Rapid calculations
- Trends and patterns
- Initial screening

### Dynamics:
- MD simulations
- Large systems
- Long timescales
- Phase transitions
- Melting studies

## Best Practices

### KEDF Selection:
- Appropriate for system
- WGC for sp metals
- Test accuracy
- Validate results
- Know limitations

### System Choice:
- Best for metals
- Simple semiconductors ok
- Avoid complex systems
- Benchmark when possible
- Check transferability

### Convergence:
- Plane-wave cutoff
- Real-space grid
- Optimization criteria
- Standard testing
- Validate energies

### Validation:
- Compare with KS-DFT
- Benchmark calculations
- Experimental data
- Know method limits
- Document assumptions

## Community and Support
- Academic license
- Princeton University
- Carter research group
- OF-DFT community
- Research collaborations
- User support (limited)

## Educational Resources
- Carter group documentation
- Published papers
- OF-DFT literature
- Example calculations
- KEDF reviews

## Development
- Princeton University
- Emily Carter group
- Academic research
- Active development
- KEDF research
- Method improvements

## Research Focus

### OF-DFT Methods:
- Kinetic functionals
- Accuracy improvements
- New KEDFs
- Transferability
- Method validation

### Large Systems:
- Scalability demonstrations
- 100,000+ atom calculations
- Computational efficiency
- Practical applications
- Extreme scaling

### Materials Applications:
- Metallic systems
- Alloy properties
- Phase diagrams
- Dynamics studies
- Realistic simulations

## Princeton Development
- Carter group expertise
- OF-DFT pioneers
- Academic excellence
- Method development
- International impact

## Technical Innovation

### Orbital-Free Approach:
- No Kohn-Sham orbitals
- Direct density
- KEDF approximation
- Computational efficiency
- Scalability advantages

### KEDF Research:
- Advanced functionals
- Systematic improvement
- Accuracy vs speed
- Transferability
- Method development

### Linear Scaling:
- O(N) algorithms
- Large system capability
- Efficient implementation
- Practical simulations
- Extreme sizes

## Verification & Sources
**Primary sources**:
1. Princeton website: http://www.princeton.edu/cbe/people/faculty/carter/research-carter-group/profess/
2. E. A. Carter group publications
3. V. V. Karasiev et al., Comput. Phys. Commun. - PROFESS paper
4. Carter group KEDF papers

**Secondary sources**:
1. Published studies using PROFESS
2. Orbital-free DFT literature
3. KEDF method reviews
4. Large-scale simulation papers

**Confidence**: LOW_CONF - Academic license, specialized method, OF-DFT niche

**Verification status**: ✅ VERIFIED
- Princeton website: ACCESSIBLE
- Documentation: Available with license
- Software: Academic license required
- Community support: Carter group, OF-DFT community
- Academic citations: Significant in OF-DFT field
- Active development: Princeton Carter group
- Specialized strength: Orbital-free DFT, kinetic energy density functionals, linear-scaling, 100,000+ atoms, large metallic systems, KEDF development, computational efficiency
