# TopMod

## Official Resources
- Homepage: https://www.lct.jussieu.fr/pagesperso/silvi/topmod_english.html
- Documentation: https://www.lct.jussieu.fr/pagesperso/fuster/TOPMOD/topmod-manual.pdf
- Source: Available from homepage
- License: Academic use

## Overview
TopMod is a FORTRAN package for topological analysis of the Electron Localization Function (ELF). It calculates the ELF on a 3D grid, assigns basins, and computes basin populations and their variances. This enables quantitative analysis of chemical bonding in terms of core, bonding, and lone pair basins.

**Scientific domain**: ELF topology, chemical bonding, electron localization
**Target user community**: Researchers studying chemical bonding through electron localization

## Theoretical Methods
- Electron Localization Function (ELF)
- Topological analysis (gradient field)
- Basin assignment and integration
- Population and variance calculation
- Synaptic order classification

## Capabilities (CRITICAL)
- **ELF Calculation**: 3D grid evaluation
- **Basin Analysis**: Core, bonding, lone pair basins
- **Population**: Basin electron populations
- **Variance**: Population fluctuation analysis
- **Synaptic Order**: Basin connectivity
- **Visualization**: Basin boundaries

**Sources**: TopMod documentation, Silvi group publications

## Key Strengths

### ELF Specialization:
- Complete ELF analysis
- Basin classification
- Population statistics
- Quantitative bonding

### Established Method:
- Well-validated
- Extensive literature
- Standard methodology
- Benchmark results

### Gaussian Compatible:
- Direct wfn input
- Standard format support
- Easy workflow
- Gaussian 92-16

## Inputs & Outputs
- **Input formats**:
  - Gaussian wfn files
  - Standard wavefunction format
  
- **Output data types**:
  - ELF values on grid
  - Basin assignments
  - Population values
  - Variance statistics

## Installation
```bash
# Download from homepage
tar -xvf topmod.tar
cd topmod
make
```

## Usage Examples
```bash
# Prepare input file
# Run TopMod
./topmod < input.dat

# Input file format:
molecule.wfn    # wavefunction file
80 80 80        # grid points
0.1             # step size
```

## Performance Characteristics
- **Speed**: Efficient grid evaluation
- **Memory**: Depends on grid size
- **Accuracy**: Well-established method

## Limitations & Known Constraints
- **FORTRAN**: Older codebase
- **Grid-based**: Resolution dependent
- **Documentation**: Limited online resources
- **Visualization**: Requires external tools

## Comparison with Other Tools
- **vs Multiwfn**: TopMod specialized ELF, Multiwfn broader
- **vs Critic2**: TopMod ELF focus, Critic2 general topology
- **vs DGrid**: Different ELF implementations
- **Unique strength**: Original ELF topology implementation

## Application Areas
- Chemical bond characterization
- Lone pair identification
- Hypervalent bonding
- Metallic bonding analysis
- Electron pair domains

## Best Practices
- Use sufficient grid resolution
- Verify basin connectivity
- Compare populations with expectations
- Cross-validate with other methods

## Community and Support
- LCT Paris group
- Academic distribution
- Published methodology
- Bernard Silvi (developer)

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.lct.jussieu.fr/pagesperso/silvi/topmod_english.html
2. S. Noury et al., Comput. Chem. 23, 597 (1999)
3. B. Silvi, A. Savin, Nature 371, 683 (1994)

**Confidence**: VERIFIED - Established academic software

**Verification status**: ✅ VERIFIED
- Homepage: ACCESSIBLE
- Documentation: AVAILABLE (PDF)
- Source code: ACADEMIC distribution
- Developer: Bernard Silvi (LCT Paris)
- Academic citations: Well-cited
- Method: ELF topology analysis
