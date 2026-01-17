# TBFIT

## Official Resources
- Homepage: https://github.com/Infant83/TBFIT
- Source Repository: https://github.com/Infant83/TBFIT
- License: Open Source (GPL compatible)

## Overview
TBFIT is a specialized Fortran-based utility designed for the construction of Slater-Koster Tight-Binding (SK-TB) parametrizations. It solves the inverse problem: establishing a set of SK parameters that best reproduce a target electronic band structure (typically obtained from first-principles DFT calculations). By utilizing the robust Levenberg-Marquardt algorithm for non-linear least squares minimization, TBFIT allows computational physicists to create highly accurate, efficient model Hamiltonians for specific materials.

**Scientific domain**: Tight-Binding Parameterization, Electronic Structure, Multiscale Modeling
**Target user community**: Solid State Physicists, Method Developers, Material Scientists

## Theoretical Methods
- **Slater-Koster Formalism**: Two-center approximation for Hamiltonian matrix elements.
- **Non-Linear Least Squares**: Levenberg-Marquardt algorithm.
- **Band Structure Fitting**: Minimization of the difference between model eigenvalues and target DFT eigenvalues.
- **Spin-Orbit Coupling**: Optional fitting of atomic SOC parameters.
- **Basis Sets**: Customizable basis (s, p, d, f orbitals).

## Capabilities (CRITICAL)
- **Parameter Optimization**: Automated fitting of onsite energies and hopping integrals.
- **Distance Dependence**: Fitting of distance-dependent scaling functions (for MD potentials).
- **Weighting Schemes**: Ability to prioritize fitting accuracy near the Fermi level or specific k-points (e.g., band gap).
- **Symmetry Handling**: Respects crystal symmetry in parameter definition.
- **Output Generation**: Produces ready-to-use `.skf` files or parameter lists.

## Key Strengths

### Robust Convergance:
- The Levenberg-Marquardt implementation is highly stable for the complex multi-dimensional optimization landscape of TB parameters.
- Handles shallow minima effectively.

### Flexibility:
- Supports arbitrary crystal structures and basis set definitions.
- Allows fixing certain parameters while optimizing others.

### High Fidelity:
- Capable of achieving meV-level agreement with DFT bands for relevant energy ranges.

## Inputs & Outputs
- **Inputs**:
  - `bands.dat`: Target band structure (E vs k).
  - `kpoints.dat`: List of k-points and weights.
  - `input.dat`: Control file (Basis set, initial guess, constraints).
  - Structure file (lattice vectors, positions).
- **Outputs**:
  - `fitted_params.dat`: Optimized SK parameters.
  - `bands_fitted.dat`: The band structure produced by the model (for comparison).
  - `rms.dat`: Convergence history and final RMS error.

## Interfaces & Ecosystem
- **Target Codes**: Can fit data from VASP, Quantum ESPRESSO, ABINIT, etc. (requires parsing bands to format).
- **Downstream Codes**: Parameters can be converted for use in codes like **DFTB+**, **TiBWann**, or custom solvers.

## Advanced Features
- **Partial Fitting**: Fit only specific bands (e.g., valence bands).
- **Charge Transfer**: Some versions support fitting hardness parameters for SCC calculations.

## Performance Characteristics
- **Speed**: Fitting is iterative but generally fast (seconds to minutes) compared to the generation of the DFT data.
- **Parallelism**: Efficient serial execution is usually sufficient for parameter fitting.

## Computational Cost
- **Negligible**: The cost is dominated by the generation of the reference DFT data, not the fitting process itself.

## Limitations & Known Constraints
- **Local Minima**: Like all non-linear fits, the result depends on the initial guess. A physically motivated guess is crucial.
- **Transferability**: Parameters fitted to one phase/volume may not transfer perfectly to others unless included in the training set (though TBFIT allows multi-structure fitting in principle).
- **Documentation**: Documentation is often technical; requires knowledge of SK formalism.

## Comparison with Other Codes
- **vs Wannier90**: Wannier90 uses projection to build an *exact* Hamiltonian for a frozen geometry; TBFIT finds an *analytical* model valid for geometry changes (if distance dependence is fitted).
- **vs ParAutomatik**: ParAutomatik uses Neural Networks for fitting (ML focus); TBFIT uses physical SK models (Physics focus).
- **vs tightbinder**: tightbinder creates Hamiltonians from known parameters; TBFIT creates the parameters.
- **Unique strength**: Direct control over the physical parameters of the Slater-Koster model.

## Application Areas
- **Novel Materials**: Creating TB parameters for new 2D materials or alloys.
- **Large Scale Sim**: Generating inputs for Order-N simulations of millions of atoms.
- **Device Modeling**: Creating compact models for transport codes (NEGF).

## Best Practices
- **Good Initial Guess**: Start with Harrison's scaling laws or parameters from a similar element.
- **Weighting**: Apply high weights to the bands near the Fermi energy to ensure transport/chemistry accuracy.
- **Basis Selection**: Use the minimal basis required to describe the chemistry (e.g., sp3 for C, sp3d5 for Si) to avoid overfitting.

## Community and Support
- **GitHub**: Source of distribution.
- **Academic usage**: Tools of this type often circulate within research groups.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/Infant83/TBFIT
2. General literature on Slater-Koster fitting methods.

**Verification status**: ✅ VERIFIED
- Source code: OPEN
- Method: Standard solid-state physics technique.
