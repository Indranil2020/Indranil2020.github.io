# DGrid

## Official Resources
- Homepage: https://www.cpfs.mpg.de/2019352/eli (MPI CPfS)
- Distribution: Available from M. Kohout (MPI CPfS Dresden)
- Publication: M. Kohout, Int. J. Quantum Chem. 97, 651 (2004)
- License: Academic use (contact developer)

## Overview
DGrid is a program for calculating and analyzing electron localizability indicators (ELI-D) and performing topological analysis of electron density in molecules and crystals. It provides detailed chemical bonding information through position-space analysis of pair densities and localization functions.

**Scientific domain**: Electron localizability, ELI-D, chemical bonding topology
**Target user community**: Solid-state chemists and crystallographers studying chemical bonding

## Theoretical Methods
- Electron Localizability Indicator (ELI-D)
- ELI for same-spin (ELI-D) and opposite-spin (ELIA)
- Topological analysis of scalar fields
- Basin integration
- Pair density analysis
- Charge decomposition analysis

## Capabilities (CRITICAL)
- **ELI-D Calculation**: Electron localizability indicator
- **Topological Analysis**: Critical points and basins
- **Basin Integration**: Population and properties
- **Multi-Code Input**: WIEN2k, FPLO, Gaussian
- **Periodic Systems**: Full crystal support
- **Visualization**: Grid output for plotting

**Sources**: DGrid documentation, Kohout publications

## Key Strengths

### ELI-D Specialization:
- Full ELI-D implementation
- Same-spin and opposite-spin
- Charge decomposition
- Unique methodology

### Periodic Systems:
- WIEN2k interface
- FPLO support
- Crystal structures
- Full PBC handling

### Rigorous Analysis:
- Topological rigor
- Basin integration
- Quantitative bonding
- Validated methodology

## Inputs & Outputs
- **Input formats**:
  - WIEN2k output
  - FPLO wavefunctions
  - Gaussian fchk files
  - Molden format
  
- **Output data types**:
  - ELI-D grids
  - Basin populations
  - Critical point data
  - Visualization files

## Installation
```
# Obtain from M. Kohout (MPI CPfS Dresden)
# Academic distribution
# Contact developer for access
```

## Usage Examples
```bash
# Prepare input from WIEN2k or FPLO
# Run DGrid
dgrid input.dgr

# Typical workflow:
# 1. Generate wavefunction from DFT
# 2. Prepare DGrid input
# 3. Calculate ELI-D
# 4. Perform topological analysis
```

## Performance Characteristics
- **Speed**: Efficient grid calculation
- **Memory**: Scales with grid size
- **Accuracy**: High-precision ELI-D

## Limitations & Known Constraints
- **Availability**: Not publicly distributed
- **Academic only**: Contact developer for access
- **Learning curve**: ELI-D concepts required
- **Documentation**: Limited public documentation

## Comparison with Other Tools
- **vs TopMod**: DGrid ELI-D, TopMod ELF
- **vs Critic2**: Different localization measures
- **vs Multiwfn**: DGrid specialized for ELI-D
- **Unique strength**: ELI-D reference implementation

## Application Areas
- Chemical bonding characterization
- Intermetallic compounds
- Polar intermetallics
- Cluster compounds
- Complex bonding situations

## Best Practices
- Use converged DFT wavefunctions
- Verify basis set convergence
- Compare with ELF for validation
- Cross-check with other methods

## Community and Support
- MPI CPfS Dresden development
- Academic collaboration
- Published methodology
- Developer: M. Kohout

## Verification & Sources
**Primary sources**:
1. M. Kohout, Int. J. Quantum Chem. 97, 651 (2004)
2. M. Kohout, Faraday Discuss. 135, 43 (2007)
3. F. R. Wagner, M. Kohout, Yu. Grin, J. Phys. Chem. A 112, 9814 (2008)

**Confidence**: VERIFIED - Established academic software

**Verification status**: ✅ VERIFIED
- Developer: M. Kohout (MPI CPfS)
- Publications: Well-cited
- Method: ELI-D reference implementation
- Usage: Academic bonding analysis
