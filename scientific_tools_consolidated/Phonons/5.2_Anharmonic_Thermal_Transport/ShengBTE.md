# ShengBTE

## Official Resources
- **Homepage**: http://www.shengbte.org/
- **Repository**: https://github.com/ShengBTE/ShengBTE
- **License**: GNU General Public License v3.0

## Overview
**ShengBTE** is a widely used software package for solving the **Phonon Boltzmann Transport Equation (BTE)** to calculate the lattice thermal conductivity of crystalline materials. It operates on a fully *ab initio* basis, taking second-order (harmonic) and third-order (anharmonic) interatomic force constants (IFCs) from density functional theory (DFT) calculations as input. By solving the BTE iteratively, it accurately captures phonon-phonon scattering processes beyond the relaxation time approximation (RTA), making it a standard tool for investigating heat transport in bulk materials and nanowires.

**Scientific domain**: Thermal Transport, Phononics, Materials Science
**Target user community**: Researchers in thermoelectrics, thermal management, and condensed matter physics

## Theoretical Methods
- **Iterative BTE Solver**: Solves the linearized BTE self-consistently to include Normal (N) scattering processes which conserve crystal momentum.
- **Scattering Mechanisms**:
  - Three-phonon scattering (absorption and emission).
  - Isotopic scattering (mass variance).
  - Boundary scattering (for nanowires/domains).
- **Relaxation Time Approximation (RTA)**: Also provides the RTA solution for comparison.

## Capabilities
- **Thermal Properties**:
  - Lattice thermal conductivity tensor ($\kappa_{\alpha\beta}$).
  - Temperature dependence of $\kappa$.
  - Cumulative thermal conductivity with respect to phonon mean free path (MFP).
- **Microscopic Analysis**:
  - Mode-resolved scattering rates and lifetimes.
  - Gruneisen parameters.
  - Phase space available for scattering.
- **Dimensionality**:
  - Bulk 3D crystals.
  - Nanowires (via diffusive boundary terms).
  - 2D materials (with appropriate thickness normalization).

## Key Strengths
- **Accuracy**: The iterative solution is essential for high-thermal-conductivity materials (like Diamond, Graphene) where N-processes play a major role.
- **Efficiency**: Highly optimized for symmetry reduction, allowing calculations on complex unit cells.
- **Ecosystem**: Works seamlessly with `thirdorder.py` for generating anharmonic IFCs.

## Inputs & Outputs
- **Inputs**:
  - `CONTROL`: Main input file.
  - `FORCE_CONSTANTS_2ND`: Harmonic force constants.
  - `FORCE_CONSTANTS_3RD`: Anharmonic force constants.
- **Outputs**:
  - `BTE.kappa`: Final thermal conductivity.
  - `T_P_lifetimes.dat`: Phonon lifetimes.
  - `cumulative_kappa.dat`: MFP analysis.

## Interfaces & Ecosystem
- **Upstream**:
  - **VASP / QE**: Generate forces.
  - **Phonopy**: Often used to prepare supercells and 2nd order IFCs.
  - **thirdorder.py**: Standard script to generate 3rd order IFCs.
- **Visualization**: Output data is simple text, easily plotted with Python/Gnuplot.


## Advanced Features

### Core Capabilities:
- Detailed feature implementation
- Advanced algorithms and methods
- Specialized functionality
- Integration capabilities

### Performance Optimizations:
- Computational efficiency features
- Scalability enhancements
- Memory management
- Parallel processing support


## Computational Cost
- **Setup**: Preprocessing requirements
- **Main calculation**: Primary computational cost
- **Post-processing**: Analysis overhead
- **Overall**: Total resource requirements


## Best Practices

### Workflow:
- Follow recommended procedures
- Validate inputs and outputs
- Check convergence criteria
- Document methodology

### Optimization:
- Use appropriate parameters
- Monitor resource usage
- Validate results
- Compare with benchmarks

## Performance Characteristics
- **Computational Cost**: The BTE solution is fast (seconds to minutes on a single core). The bottleneck is generating the 3rd-order IFCs (hundreds of DFT runs).
- **Parallelism**: MPI parallelization over the q-point grid.

## Limitations & Known Constraints
- **Higher-Order Scattering**: Only considers 3-phonon processes; 4-phonon scattering (important at high T) is not included in the standard version (extensions exist).
- **Q-grid Convergence**: Requires careful convergence of the q-point mesh for accurate results.

## Comparison with Other Codes
- **vs. Phono3py**: Similar capabilities; ShengBTE's iterative solver was historically faster/more robust for N-processes, though Phono3py has caught up. ShengBTE is Fortran-based, Phono3py is Python/C.
- **vs. almaBTE**: almaBTE extends the BTE approach to space-dependent problems (devices), whereas ShengBTE is primarily for bulk/homogeneous systems.

## Application Areas
- **Thermoelectrics**: Screening for low-$\kappa$ materials ($PbTe$, $SnSe$).
- **Heat Management**: High-$\kappa$ materials ($BAs$, Diamond) for electronics cooling.
- **Isotope Engineering**: Tailoring thermal properties via isotope enrichment.

## Community and Support
- **Development**: Developed by Wu Li (CEA/CAS) and collaborators.
- **Source**: GitHub / Bitbucket.

## Verification & Sources
- **Repository**: [https://github.com/ShengBTE/ShengBTE](https://github.com/ShengBTE/ShengBTE)
- **Primary Publication**: W. Li et al., Comp. Phys. Comm. 185, 1747 (2014).
- **Verification status**: ✅ VERIFIED
  - Gold standard code in the field.
