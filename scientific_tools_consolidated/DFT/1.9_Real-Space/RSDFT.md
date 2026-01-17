# RSDFT (Real-Space Density Functional Theory)

## Official Resources
- Homepage: http://rsdft.jp/
- Source Repository: https://github.com/j-iwata/RSDFT
- License: Open Source (GPLv3)

## Overview
RSDFT is a high-performance ab initio density functional theory code based on the real-space finite-difference method. It is explicitly designed for massively parallel computing architectures, enabling first-principles calculations on systems containing over 100,000 atoms. The code avoids the global communication bottlenecks of Fast Fourier Transforms (FFTs) typical in plane-wave codes, making it highly scalable on supercomputers like the K computer and Fugaku.

**Scientific domain**: Large-scale materials science, Nanostructures, Surfaces, Interfaces
**Target user community**: HPC researchers, large-scale DFT practitioners

## Theoretical Methods
- Real-Space Finite-Difference Method
- Kohn-Sham Density Functional Theory
- Pseudopotentials (Norm-conserving and Ultrasoft)
- Real-space grid discretization
- High-order finite difference stencils (e.g., higher-order central differences)
- Conjugate Gradient and RMM-DIIS eigensolvers
- Projector Augmented Wave (PAW) method (experimental/supported in versions)

## Capabilities
- **System Size**: Routine calculations for 10,000 to 100,000+ atoms.
- **Parallelism**: Scalable to >80,000 nodes (Gordon Bell Prize Winner 2011).
- **Electronic Structure**: Ground state energy, forces, stresses.
- **Molecular Dynamics**: Born-Oppenheimer MD for large systems.
- **Boundary Conditions**: Flexible (Isolated, Periodic, Wire, Slab).

## Key Strengths

### Massively Parallel Scalability:
- No FFT bottleneck.
- Domain decomposition allows linear scaling with respect to processors for fixed system size per core.
- Demonstrated on exascale-class hardware.

### Large-Scale Calculations:
- Can handle silicon nanowires with 100k+ atoms.
- Ideal for complex nanostructures where periodic boundary conditions of plane-waves are artificial.

### Flexible Boundary Conditions:
- Naturally handles non-periodic systems without supercell approximation artifacts (charged systems, dipoles).

## Inputs & Outputs
- **Inputs**:
  - Grid spacing
  - Finite difference order
  - Atomic coordinates
  - Pseudopotential files
- **Outputs**:
  - Total energies and forces
  - Charge density distribution (cube files)
  - Wavefunctions (real-space grid)

## Interfaces & Ecosystem
- **Python**: Prototyping environment available in Python.
- **HPC Integration**: Optimized for MPI/OpenMP hybrid parallelism.
- **MateriApps**: Listed and supported via the MateriApps ecosystem.

## Advanced Features
- **GPU Acceleration**: Ports available for GPU-accelerated clusters.
- **Order-N capability**: Methods for linear scaling electronic structure (Chebyshev filtering).

## Performance Characteristics
- **Speed**: Superior to Plane-Wave codes for very large systems (>1000 atoms) on massive core counts.
- **Accuracy**: Systematic convergence with grid spacing (h) and stencil order.
- **Efficiency**: Excellent weak scaling.

## Computational Cost
- **High**: For small systems (under 100 atoms), overhead is higher than VASP/QE.
- **Low**: For large systems (>5000 atoms), significantly cheaper/faster than PW codes.

## Limitations & Known Constraints
- **Pseudopotentials**: requires real-space optimized potentials for best efficiency.
- **Grid anisotropy**: \"Egg-box\" effect possible if grid is too coarse (breaking translation symmetry).
- **Documentation**: English documentation can be sparse compared to VASP; website is partly in Japanese.

## Comparison with Other Codes
- **vs PARSEC**: Both are real-space; RSDFT focuses more on massive HPC scalability.
- **vs ONETEP**: ONETEP attempts O(N) via Wannier functions; RSDFT uses domain decomposition of the grid.
- **vs Plane-Wave**: RSDFT eliminates FFT communication, winning at extreme scale.
- **Unique strength**: Gordon Bell Prize-winning scalability for 100k atom systems.

## Application Areas
- **Nanowires**: Electronic properties of realistic diameter wires.
- **Biomolecules**: Large proteins in solvent (implicit/explicit).
- **Defects**: Dilute defects in very large supercells.

## Best Practices
- **Grid Convergence**: Test h-grid spacing carefully to avoid egg-box errors.
- **Parallel Layout**: Match domain decomposition to physical system shape.
- **Pseudopotentials**: Use soft potentials where possible to allow coarser grids.

## Community and Support
- **MateriApps**: Integration with the Japanese materials science setup.
- **Development Team**: Based at U-Tokyo and RIKEN.

## Verification & Sources
**Primary sources**:
1. Official Website: http://rsdft.jp/
2. GitHub: https://github.com/j-iwata/RSDFT
3. J. Hasegawa et al., "First-Principles Calculations of 100,000-Atom Silicon Nanowires", SC11 (Gordon Bell Prize).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GPLv3)
- Performance: Confirmed high-impact HPC code.
