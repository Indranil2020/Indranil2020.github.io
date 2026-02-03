# Phono3py

## Official Resources
- **Homepage**: https://phonopy.github.io/phono3py/
- **Repository**: https://github.com/phonopy/phono3py
- **License**: BSD-3-Clause

## Overview
**Phono3py** is a premier open-source software package for calculating **phonon-phonon interactions** and **lattice thermal conductivity** from first principles. Built upon the foundation of **Phonopy**, it enables the computation of third-order force constants (anharmonicity) and solves the Boltzmann Transport Equation (BTE) to predict thermal properties of crystals. It is widely used in the materials science community due to its robustness, ease of use, and tight integration with major DFT codes.

**Scientific domain**: Lattice Dynamics, Thermal Transport, Anharmonicity
**Target user community**: Materials scientists studying thermal properties and phonon lifetimes

## Theoretical Methods
- **Finite Displacement Method**: Constructs third-order force constants (FC3) by displacing atoms in a supercell.
- **Phonon BTE**: Solves the linearized BTE using the Relaxation Time Approximation (RTA) or the full iterative self-consistent method (LBTE).
- **Phonon Self-Energy**: Calculates the imaginary part of the self-energy (linewidth) due to 3-phonon scattering.
- **Joint Density of States (JDOS)**: Analyzes the phase space available for scattering processes.

## Capabilities
- **Thermal Transport**:
  - Lattice thermal conductivity tensor ($\kappa_L$).
  - Accumulation of $\kappa_L$ with respect to frequency and mean free path.
- **Anharmonic Properties**:
  - Phonon lifetimes ($\tau$) and linewidths ($\Gamma$).
  - Gruneisen parameters (from FC3).
  - Frequency shifts due to thermal expansion (quasi-harmonic approx via Phonopy) or self-energy (real part).
- **Wigner Transport**: Recent implementation of the Wigner transport formalism for handling complex crystals with particle-like and wave-like heat transport.

## Key Strengths
- **Integration**: Works directly with **Phonopy** inputs and outputs, providing a unified workflow for harmonic and anharmonic properties.
- **Documentation**: Extensive and high-quality documentation with tutorials.
- **Interface**: Python API ("phonopy-like") allows for flexible scripting and post-processing; C-extension for performance.
- **Algorithmic Depth**: Supports tetrahedron method for accurate Brillouin zone integration.

## Inputs & Outputs
- **Inputs**:
  - Structure file (POSCAR, unitcell).
  - Force sets from supercell DFT calculations (very large number of runs).
  - `disp_fc3.yaml`: Displacement dataset.
- **Outputs**:
  - `kappa-m*.hdf5`: Detailed thermal conductivity data.
  - `fc3.hdf5`: Third-order force constants.
  - `gammas`: Phonon linewidths.

## Interfaces & Ecosystem
- **DFT Codes**: VASP, Quantum ESPRESSO, CRYSTAL, TURBOMOLE, CASTEP, Wien2k, etc.
- **Phonopy**: Core dependency.
- **HDF5**: Heavy data (FC3) stored in HDF5 format for efficiency.


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
- **Computational Cost**: The bottleneck is the generation of forces for supercells with 2-atom displacements (scaling roughly as $N_{atoms}^2$ or $N_{atoms}^3$ depending on cutoff).
- **Solver Speed**: The BTE solver is OpenMP parallelized and efficient for moderate k-grids.
- **Memory**: Storing the full dense FC3 matrix can be memory intensive for large unit cells.

## Limitations & Known Constraints
- **Cost**: Requires significantly more computational resources than harmonic Phonopy calculations.
- **4-Phonon**: Standard version includes only 3-phonon processes (4-phonon extensions are experimental or separate).

## Comparison with Other Codes
- **vs. ShengBTE**: ShengBTE uses 3rd order constants too (often from `thirdorder.py`); Phono3py has its own displacement generator. Phono3py is more "Pythonic" and integrated with Phonopy's workflow.
- **vs. almaBTE**: almaBTE focuses on BTE solvers for devices/multiscale; Phono3py is primarily for bulk properties.

## Application Areas
- **Thermoelectric Materials**: Prediction of $zT$ limits.
- **Thermal Barrier Coatings**: Search for ultra-low conductivity oxides.
- **Fundamental Physics**: Study of phonon hydrodynamics and Wigner transport regimes.

## Community and Support
- **Development**: Approached and maintained by Atsushi Togo (NIMS, Japan).
- **Forum**: Active user mailing list.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/phonopy/phono3py](https://github.com/phonopy/phono3py)
- **Primary Publication**: A. Togo et al., Phys. Rev. B 91, 094306 (2015).
- **Verification status**: ✅ VERIFIED
  - Widely used and trusted community code.
