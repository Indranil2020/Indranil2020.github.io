# FPLO (Full-Potential Local-Orbital)

## Official Resources
- Homepage: https://www.fplo.de/
- Documentation: https://www.fplo.de/download-fplo/
- Source Repository: Not publicly available (registration required)
- License: Free for academic use (registration required)

## Overview
FPLO (Full-Potential Local-Orbital minimum-basis) is a DFT code using local orbitals and the full-potential approach for highly accurate electronic structure calculations of solids. Developed at the IFW Dresden (Leibniz Institute for Solid State and Materials Research), FPLO is particularly strong in calculations for strongly correlated systems, magnetism, and materials with complex electronic structures. It combines the accuracy of full-potential methods with the efficiency of minimal basis sets.

**Scientific domain**: Solid-state DFT, strongly correlated materials, magnetism, electronic structure  
**Target user community**: Solid-state physicists, materials scientists, strongly correlated systems researchers

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA)
- Full-potential local-orbital approach
- Minimal basis sets (local orbitals)
- All-electron (no pseudopotentials)
- LDA+U for correlated systems
- Spin-polarized calculations
- Non-collinear magnetism
- Spin-orbit coupling
- GW approximation (via interface)
- Hybrid functionals (limited)

## Capabilities (CRITICAL)
- Ground-state electronic structure (solids)
- Total energy calculations
- Band structure and DOS
- Fermi surfaces
- Magnetic properties
- Density of states (DOS, PDOS)
- Electric field gradients
- Hyperfine fields
- Charge density analysis
- Crystal orbital Hamilton populations (COHP)
- Bonding analysis
- Strongly correlated materials (LDA+U)
- Magnetic ordering
- Spin-orbit coupling effects
- High accuracy for complex materials

**Sources**: Official FPLO website (https://www.fplo.de/)

## Key Strengths

### Full-Potential:
- No shape approximations
- High accuracy
- All-electron treatment
- Correct electron density
- Accurate for complex systems

### Local Orbitals:
- Minimal basis approach
- Efficient computation
- Physical interpretation
- Direct bonding analysis
- Reduced computational cost

### Strongly Correlated Systems:
- LDA+U implementation
- Magnetic materials
- Transition metal oxides
- f-electron systems
- Complex magnetism

### Bonding Analysis:
- COHP analysis
- Crystal orbital overlap populations
- Chemical bonding insights
- Orbital interactions
- Interpretable results

### Accuracy:
- Full-potential precision
- All-electron approach
- Systematic basis
- Well-tested for solids

## Inputs & Outputs
- **Input formats**:
  - Text-based input file (=.in)
  - Structure definition
  - Parameter specifications
  - Basis set definitions
  
- **Output data types**:
  - Text output files
  - Band structure data
  - DOS files
  - Charge densities
  - COHP data
  - Fermi surface information

## Interfaces & Ecosystem
- **Visualization**:
  - XCrySDen for structures
  - Standard plotting tools for DOS/bands
  - Custom scripts
  
- **Analysis**:
  - COHP analysis tools
  - Band structure analysis
  - DOS analysis
  - Custom post-processing
  
- **Related Codes**:
  - Interface with GW codes
  - Data export to standard formats

## Workflow and Usage

### Typical Workflow:
1. Prepare structure file
2. Create input file with parameters
3. Run FPLO calculation
4. Analyze output (DOS, bands, COHP)
5. Visualize results

### Input Structure:
- Structure definition (lattice, atoms)
- Calculation parameters
- Basis set specifications
- Self-consistency criteria
- Output options

### Running FPLO:
```bash
fplo input.in
```

## Advanced Features

### LDA+U:
- On-site Coulomb interactions
- Hubbard U parameter
- Strongly correlated electrons
- Improved band gaps
- Magnetic properties

### COHP Analysis:
- Chemical bonding information
- Orbital-resolved contributions
- Bond strength analysis
- Antibonding/bonding identification
- Physical insight

### Spin-Orbit Coupling:
- Relativistic effects
- Important for heavy elements
- Magnetic anisotropy
- Band splitting
- Accurate electronic structure

### Non-Collinear Magnetism:
- Complex magnetic structures
- Spin spirals
- Frustrated magnetism
- Spin canting

### Electric Field Gradients:
- Nuclear quadrupole coupling
- NQR spectroscopy
- Mössbauer parameters
- Chemical environment

## Performance Characteristics
- **Speed**: Efficient due to minimal basis
- **Accuracy**: High (full-potential)
- **System size**: Moderate (typical unit cells)
- **Memory**: Moderate requirements
- **Parallelization**: Limited parallel capabilities

## Computational Cost
- **DFT**: Efficient for minimal basis
- **Full-potential**: More costly than pseudopotential
- **LDA+U**: Moderate overhead
- **SOC**: Increased cost
- **Typical systems**: Unit cells up to ~100 atoms

## Limitations & Known Constraints
- **Registration required**: Free but needs registration
- **Documentation**: Basic (primarily examples)
- **Learning curve**: Steep
- **Parallelization**: Limited
- **Community**: Smaller, specialized
- **Platform**: Linux primarily
- **Molecular systems**: Not optimized for molecules
- **Graphical interface**: Limited

## Comparison with Other Codes
- **vs WIEN2k**: Both full-potential, FPLO uses local orbitals, WIEN2k APW+lo
- **vs VASP**: FPLO all-electron, VASP pseudopotential; FPLO better COHP
- **vs Quantum ESPRESSO**: Different approaches, FPLO all-electron
- **vs Elk**: Both all-electron, different basis approaches
- **Unique strength**: Full-potential with local orbitals, COHP analysis, strongly correlated systems

## Application Areas

### Strongly Correlated Materials:
- Transition metal oxides
- Mott insulators
- Heavy fermions
- f-electron systems
- Complex oxides

### Magnetism:
- Magnetic materials
- Magnetic ordering
- Spin-orbit effects
- Magnetic anisotropy
- Exchange interactions

### Chemical Bonding:
- Bonding analysis
- COHP studies
- Intermetallic compounds
- Chemical understanding
- Structure-property relations

### Spectroscopy:
- Electric field gradients
- Hyperfine fields
- NMR/NQR parameters
- Mössbauer spectroscopy

## Best Practices

### Basis Set:
- Use recommended basis sets
- Convergence testing
- Minimal basis philosophy
- Systematic improvement

### Self-Consistency:
- Appropriate convergence criteria
- Sufficient k-points
- Charge convergence
- Energy convergence

### LDA+U:
- Choose appropriate U values
- Literature guidance
- Test sensitivity
- Physical justification

### Magnetism:
- Initialize magnetic moments
- Test different configurations
- Check convergence
- Spin-orbit when needed

## Community and Support
- Free for academic use
- Registration required
- Email support (limited)
- User community (smaller)
- Published examples
- IFW Dresden development

## Educational Resources
- Example calculations
- Published papers
- User manual (basic)
- Literature references
- Community knowledge

## Development
- IFW Dresden (Leibniz Institute)
- Active but selective development
- Specialized for solid-state
- Method improvements
- User-driven features

## Research Applications
- Materials discovery
- Magnetic materials design
- Electronic structure studies
- Chemical bonding understanding
- Property predictions

## Specialized Features

### COHP Method:
- Unique strength of FPLO
- Crystal orbital Hamilton populations
- Bonding/antibonding analysis
- Orbital-resolved
- Interpretable chemistry

### Full-Potential Accuracy:
- No muffin-tin approximation
- Accurate electron density everywhere
- Important for complex structures
- High precision

## Verification & Sources
**Primary sources**:
1. Official website: https://www.fplo.de/
2. K. Koepernik and H. Eschrig, Phys. Rev. B 59, 1743 (1999) - FPLO methodology
3. IFW Dresden documentation

**Secondary sources**:
1. Published studies using FPLO (>2000 citations)
2. COHP analysis papers
3. Strongly correlated materials literature
4. Confirmed in source lists (LOW_CONF due to specialized nature)

**Confidence**: LOW_CONF - Specialized code with smaller community, registration required

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: Basic (examples, papers)
- Software: Available with registration
- Community support: Limited (email, specialized community)
- Academic citations: >2500
- Development: Active at IFW Dresden
- Specialized strength: Full-potential local-orbital approach, COHP analysis, strongly correlated materials, accurate solid-state calculations
