# DDEC (Density Derived Electrostatic and Chemical)

## Official Resources
- Homepage: https://sourceforge.net/projects/ddec/
- Documentation: https://sourceforge.net/projects/ddec/files/
- Source Repository: https://sourceforge.net/projects/ddec/
- License: MIT License

## Overview
DDEC (Density Derived Electrostatic and Chemical) is a method and software code (`chargemol`) for computing net atomic charges, atomic spin moments, and bond orders from quantum mechanical charge density distributions. The DDEC6 method is designed to reproduce the electrostatic potential outside the molecular charge distribution while maintaining chemical transferability and spherical averaging convergence.

**Scientific domain**: Charge analysis, bond order analysis, electrostatics  
**Target user community**: Computational chemists, materials scientists

## Theoretical Methods
- Density Derived Electrostatic and Chemical (DDEC) charge partitioning
- DDEC6 method (latest version)
- Atomic spin moment calculation
- Effective Bond Order (EBO)
- Overlap populations
- Electrostatic potential fitting

## Capabilities (CRITICAL)
- Calculation of net atomic charges (NACs) from electron density
- Calculation of atomic spin moments (ASMs)
- Computation of bond orders for all atom pairs
- Handling of non-periodic and periodic systems
- Reproduces electrostatic potential accurately
- Stable against basis set variations
- Optimized for speed (OpenMP parallelization)

**Sources**: DDEC documentation, RSC Adv. 6, 27724 (2016)

## Key Strengths

### DDEC6 Method:
- Reproduces electrostatic potential
- Chemically transferable charges
- Robust for charged systems
- All elements supported

### Bond Order Calculation:
- Comprehensive bond orders
- Connectivity determination
- Spin moment partitioning
- Overlap populations

### Performance:
- OpenMP parallelized
- Efficient memory usage
- Large system capable
- Fast execution

## Inputs & Outputs
- **Input formats**: VASP (AECCAR0/2, CHGCAR, POTCAR), Gaussian (.wfn/.fchk), Q-Chem, etc.
- **Output data types**: `net_atomic_charges.xyz`, `overlap_populations.xyz`, `atomic_spin_moments.xyz`, statistics logs

## Interfaces & Ecosystem
- **VASP**: Direct support via charge density files
- **Gaussian/Q-Chem**: Via formatted checkpoint or cube files
- **CP2K/ONETEP**: Supported formats
- **Python**: Tools available for parsing output

## Workflow and Usage
1. Perform QM calculation (DFT).
2. Save charge density (and core density for VASP).
3. Prepare `job_control.txt` for chargemol.
4. Run `chargemol_FORTRAN_09_26_2017` (or newer executable).
5. Analyze `.xyz` output files.

## Performance Characteristics
- Highly efficient (seconds to minutes for typical systems)
- Parallelized with OpenMP
- Memory efficient

## Limitations & Known Constraints
- **Atomic densities required**: Needs reference density files
- **Input format**: Specific file format requirements
- **Documentation**: Could be more extensive
- **Learning curve**: Setup requires careful attention

## Comparison with Other Tools
- **vs Bader**: DDEC stockholder-based, Bader topology-based
- **vs Hirshfeld**: DDEC6 more robust for charged systems
- **vs Mulliken**: DDEC basis-set independent
- **Unique strength**: Combined charges and bond orders

## Application Areas
- Force field parameterization (partial charges)
- MOF/COF adsorption studies
- Reactivity prediction
- Magnetic moment distribution
- Bond characterization

## Best Practices
- Include core density for VASP (AECCAR files)
- Verify net charge matches expected
- Check bond order reasonableness
- Use appropriate atomic reference densities

## Community and Support
- Developed by Thomas Manz group (New Mexico State University)
- Hosted on SourceForge
- Citation: T. A. Manz and N. Gabaldon Limas, RSC Adv. 6, 27724 (2016)

## Verification & Sources
**Primary sources**:
1. Homepage: https://sourceforge.net/projects/ddec/
2. Publication: T. A. Manz and N. Gabaldon Limas, RSC Adv. 6, 27724 (2016)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (SourceForge)
- Documentation: AVAILABLE
- Source: OPEN (MIT)
- Development: ACTIVE (Manz Group)
- Applications: DDEC6 charges, bond orders, spin moments
