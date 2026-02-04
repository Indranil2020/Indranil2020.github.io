# Chargemol

## Official Resources
- Homepage: https://sourceforge.net/projects/ddec/
- GitHub Mirror: https://github.com/berquist/chargemol
- Documentation: Instructions PDF in download
- Publication: T. A. Manz, RSC Adv. 7, 45552 (2017)
- License: Check with authors

## Overview
Chargemol is a program for computing DDEC (Density Derived Electrostatic and Chemical) atomic charges, atomic spin moments, and bond orders from electron and spin density distributions. It implements DDEC6, the latest version of the DDEC method, providing chemically meaningful atomic charges for molecular simulations.

**Scientific domain**: Charge partitioning, bond order analysis, force field development
**Target user community**: Force field developers, molecular dynamics researchers, charge analysis

## Theoretical Methods
- DDEC6 atomic population analysis
- Iterative stockholder partitioning
- Bond order calculation (Mayer-like)
- Spin moment partitioning
- Overlap population analysis

## Capabilities (CRITICAL)
- **DDEC6 Charges**: State-of-the-art charge partitioning
- **Bond Orders**: Comprehensive bond order calculation
- **Spin Moments**: Atomic spin magnetization
- **Multiple Formats**: VASP, Gaussian, QE support
- **Periodic Systems**: Full periodic boundary support
- **Large Systems**: Efficient for thousands of atoms

**Sources**: Chargemol documentation, RSC Advances publications

## Key Strengths

### DDEC6 Method:
- Chemically meaningful charges
- Transferable across environments
- Robust for all elements
- Validated extensively

### Force Field Development:
- Ideal for parameterization
- Consistent methodology
- Spin-polarized systems
- Bond order for connectivity

### Broad Compatibility:
- VASP CHGCAR/AECCAR
- Gaussian cube files
- Quantum ESPRESSO
- Other cube formats

## Inputs & Outputs
- **Input formats**:
  - VASP CHGCAR, AECCAR0, AECCAR2
  - Gaussian cube files
  - Quantum ESPRESSO pp.x output
  
- **Output data types**:
  - DDEC6 atomic charges
  - Bond orders
  - Atomic spin moments
  - Overlap populations

## Installation
```bash
# Download from SourceForge
# Extract and compile
cd chargemol_09_26_2017
make
```

## Usage Examples
```bash
# Create job_control.txt
<net charge>
0.0
</net charge>
<periodicity along A, B, and C vectors>
.true.
.true.
.true.
</periodicity along A, B, and C vectors>
<atomic densities directory complete path>
/path/to/atomic_densities/
</atomic densities directory complete path>

# Run Chargemol
./Chargemol_09_26_2017_linux_parallel
```

## Performance Characteristics
- **Speed**: Efficient parallel implementation
- **Memory**: Scales with system size
- **Parallelization**: OpenMP support

## Limitations & Known Constraints
- **Atomic densities**: Requires reference atomic density files
- **Input format**: Specific format requirements
- **Documentation**: Could be more extensive
- **Compilation**: May require adjustments

## Comparison with Other Tools
- **vs Bader**: Different partitioning philosophy
- **vs DDEC (SourceForge)**: Chargemol is the implementation
- **vs Hirshfeld**: DDEC6 more robust for charged systems
- **Unique strength**: DDEC6 method, bond orders included

## Application Areas
- Force field parameterization
- Molecular dynamics charges
- Reactivity analysis
- Electrostatic potential fitting
- Materials property prediction

## Best Practices
- Include core density (AECCAR files for VASP)
- Verify charge neutrality
- Check bond orders against chemical intuition
- Use appropriate reference densities

## Community and Support
- SourceForge project
- Published methodology
- Developer: Thomas Manz
- Academic support

## Verification & Sources
**Primary sources**:
1. SourceForge: https://sourceforge.net/projects/ddec/
2. T. A. Manz, RSC Adv. 7, 45552 (2017) - DDEC6 bond orders
3. T. A. Manz, N. Gabaldon Limas, RSC Adv. 6, 47771 (2016) - DDEC6 charges

**Confidence**: VERIFIED - Published in RSC Advances

**Verification status**: ✅ VERIFIED
- SourceForge project: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: AVAILABLE
- Academic citations: Well-cited
- Active development: Maintained
