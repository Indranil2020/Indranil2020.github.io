# DoNOF.jl

## Official Resources
- Homepage: https://github.com/felipelewyee/DoNOF.jl
- Documentation: https://github.com/felipelewyee/DoNOF.jl/blob/main/README.md
- Source Repository: https://github.com/felipelewyee/DoNOF.jl
- License: MIT License

## Overview
DoNOF.jl is a Julia implementation of Natural Orbital Functional (NOF) theory for electronic structure calculations. It provides an alternative to traditional DFT and wavefunction methods by describing the electronic correlation through functionals of the one-particle reduced density matrix. DoNOF.jl enables calculations using PNOF functionals (PNOF5, PNOF7, GNOF) which capture static and dynamic correlation effects.

**Scientific domain**: Molecules, strongly correlated systems, multi-reference problems  
**Target user community**: Researchers studying strongly correlated systems, those needing multi-reference quality without wavefunction complexity

## Theoretical Methods
- Natural Orbital Functional (NOF) theory
- PNOF5 functional (intrapair correlation)
- PNOF7 functional (static + dynamic correlation)
- GNOF (global NOF, balanced description)
- Gaussian Type Orbital (GTO) basis sets
- Orbital optimization algorithms
- Occupation number optimization
- One-particle reduced density matrix formalism

## Capabilities (CRITICAL)
- Ground-state electronic structure via NOF
- PNOF5 calculations (intrapair electron correlation)
- PNOF7 calculations (interpair correlation)
- GNOF calculations (global NOF)
- Orbital and occupation optimization
- Total energy calculations
- Natural orbital analysis
- Fractional occupation numbers
- Strong correlation treatment
- Multi-reference character description

**Sources**: GitHub repository, M. Piris group publications

## Key Strengths

### Natural Orbital Functional Theory:
- Direct correlation functional
- Avoids wavefunction complexity
- Fractional occupations natural
- Size-consistent formulation
- Captures static correlation

### Strong Correlation Capability:
- Multi-reference problems accessible
- Bond breaking description
- Biradical systems
- Transition metal complexes
- No active space selection needed

### Julia Implementation:
- High-performance dynamic language
- Just-in-time compilation
- Interactive development
- Modern language features
- Easy to extend

### PNOF Family:
- PNOF5 for dynamic correlation
- PNOF7 adds static correlation
- GNOF balanced treatment
- Systematic improvability

## Inputs & Outputs
- **Input formats**:
  - Julia API
  - Molecular geometry specification
  - Basis set selection
  
- **Output data types**:
  - Total NOF energies
  - Natural orbital occupations
  - Natural orbitals
  - One-particle density matrix

## Interfaces & Ecosystem
- **Julia ecosystem**:
  - Julia package manager installation
  - REPL interactive use
  - Jupyter/Pluto notebook support
  
- **Integral evaluation**:
  - GaussianBasis.jl integration
  - Standard basis sets
  - Custom basis support

## Advanced Features

### NOF Variants:
- PNOF5: Intrapair correlation only
- PNOF7: Static + dynamic correlation
- GNOF: Global balanced functional
- Extensible to new functionals

### Optimization:
- Orbital coefficient optimization
- Occupation number constraints
- Convergence acceleration
- Multiple algorithms

### Analysis:
- Natural orbital populations
- Correlation indices
- Entanglement measures
- Strong correlation diagnostics

### Multi-Reference Character:
- Automatic handling
- No active space needed
- Smooth potential surfaces
- Correct dissociation limits

## Performance Characteristics
- **Speed**: Julia JIT compiled efficiency
- **Accuracy**: Multi-reference quality
- **System size**: Small to medium molecules
- **Memory**: Moderate requirements
- **Parallelization**: Julia threading

## Computational Cost
- **NOF**: Similar to coupled cluster scaling
- **Orbital optimization**: Iterative refinement
- **Convergence**: Usually 20-50 iterations
- **Typical**: Minutes to hours for molecules

## Limitations & Known Constraints
- **System size**: Small to medium molecules
- **Maturity**: Developing project
- **Periodicity**: Molecular only
- **Properties**: Energy-focused currently
- **Documentation**: Developing
- **User base**: Specialized NOF community

## Comparison with Other Codes
- **vs Traditional DFT**: DoNOF captures static correlation
- **vs CASSCF**: No active space selection needed
- **vs Coupled Cluster**: Different correlation treatment
- **vs Fermi.jl**: DoNOF specializes in NOF, Fermi.jl DFT/CC
- **Unique strength**: NOF methodology in Julia, strong correlation

## Application Areas

### Strongly Correlated Systems:
- Transition metal complexes
- Open-shell systems
- Biradicals
- Stretched bonds

### Bond Breaking:
- Dissociation curves
- Correct asymptotes
- Smooth potential surfaces
- Reaction coordinates

### Multi-Reference Problems:
- Near-degeneracy
- Avoided crossings
- Conical intersections (approximate)
- Excited state character

### Method Development:
- New NOF functional testing
- Correlation analysis
- Benchmark calculations
- Algorithm development

## Best Practices

### Getting Started:
- Install via Julia package manager
- Start with small test molecules
- Compare with reference data
- Understand NOF theory basics

### Functional Selection:
- PNOF5 for simple systems
- PNOF7 for static correlation
- GNOF for balanced description
- Test on benchmarks

### Convergence:
- Monitor energy and occupations
- Check orbital stability
- Adjust optimization parameters
- Verify convergence criteria

### Validation:
- Compare with CASPT2 for strong correlation
- Check single-reference diagnostics
- Validate on known benchmarks

## Community and Support
- Open source MIT license
- GitHub repository
- M. Piris group (San Sebastian)
- NOF method publications
- Growing Julia chemistry community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/felipelewyee/DoNOF.jl
2. M. Piris et al., NOF theory publications
3. Julia package registry

**Secondary sources**:
1. Natural orbital functional theory literature
2. Strong correlation methodology papers
3. Julia computational chemistry community

**Confidence**: VERIFIED - Active GitHub, published methodology

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Package: Julia registry
- Academic basis: M. Piris NOF methodology
- Active development: Recent commits
- Specialty: Natural Orbital Functional theory, strong correlation, Julia
