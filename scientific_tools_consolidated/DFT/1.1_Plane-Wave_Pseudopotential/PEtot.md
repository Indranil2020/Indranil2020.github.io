# PEtot

## Official Resources
- Homepage: https://psi-k.net/codes/petot (Archive/Reference)
- Documentation: Available within distribution packages (often manual.pdf)
- Source Repository: No central modern repo; historically distributed by L.-W. Wang's group (LBNL)
- License: Open Source (specifics vary by distribution, often academic/research use)

## Overview
PEtot is a plane-wave pseudopotential Density Functional Theory (DFT) code specifically designed for large-scale materials simulations. It employs norm-conserving and ultrasoft pseudopotentials and is renowned for its efficient parallelization capabilities, allowing it to scale to thousands of processors. It serves as a foundational engine for other codes, such as PWtransport for quantum transport calculations.

**Scientific domain**: Materials science, condensed matter physics, quantum transport (as backend)
**Target user community**: Researchers studying large systems, nanostructures, and needing a highly parallelizable plane-wave code

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets
- Norm-conserving pseudopotentials
- Ultrasoft pseudopotentials
- LDA, GGA functionals
- Iterative diagonalization (Conjugate Gradient, DIIS)
- Empirical pseudopotential method (EPM) capabilities (in some variants)
- Real-space multigrid options (in related specialized versions)

## Capabilities
- Ground-state electronic structure
- Large-scale system simulation (thousands of atoms)
- Geometry optimization (atomic relaxation)
- Ab initio Molecular Dynamics (AIMD)
- Band structure calculation
- Total energy and force calculations
- Massive parallelization (MPI)
- Support for isolated and periodic systems
- Electronic density generation

## Key Strengths

### Scalability and Parallelism:
- Designed for High-Performance Computing (HPC)
- Efficient MPI parallelization over G-vectors, bands, and k-points
- Scales well to thousands of cores
- Suitable for large supercomputers

### Application Areas
- **Semiconductor Physics**: Large supercell calculations for defects and dopants.
- **Nanostructures**: Quantum dots, nanowires, and isolated clusters.
- **Quantum Transport**: Electronic structure generation for NEGF codes like PWtransport.
- **Alloys**: Disorder modeling using large supercells.

### Large System Handling:
- Optimized for systems with high atom counts
- Efficient memory management
- Fast iterative solvers for large Hamiltonians

### Algorithm Versatility:
- Multiple wave function solution methods (Band-by-band, All-band CG, All-band DIIS)
- Flexibility in handling different system types (insulators, metals)

## Inputs & Outputs
- **Input formats**:
  - Main input file (control parameters)
  - Atom configuration files (coordinates)
  - Pseudopotential files (UPF or native formats)
  
- **Output data types**:
  - Standard output (energies, forces, convergence info)
  - Charge density files
  - Wavefunction files
  - Band structure data
  - Relaxed coordinates

## Interfaces & Ecosystem
- **Integrations**:
  - Foundational engine for PWtransport
  - Used in conjunction with various post-processing tools in the L.-W. Wang group ecosystem
  - Interfaces with visualization tools capable of reading standard charge/density formats

## Performance Characteristics
- **Speed**: High performance on parallel architectures due to optimized MPI communication.
- **System size**: Capable of handling thousands of atoms, making it competitive with other flagship codes for large-scale problems.
- **Parallelization**: Multi-level parallelization ensures efficient resource utilization.

## Computational Cost
- **Scaling**: Generally $O(N^3)$ with system size, but efficient prefactors due to optimized FFTS and linear algebra.
- **Memory**: Distributed memory model allows handling systems that would exceed single-node RAM.
- **Efficiency**: Plane-wave basis requires large grids for "hard" pseudopotentials, but ultrasoft/norm-conserving potentials mitigate this.

## Best Practices

### Parallelization Strategy:
- **Hybrid MPI**: Use MPI for inter-node communication; pure MPI is often standard for PEtot.
- **K-point Distribution**: For smaller systems, distribute k-points first for near-linear scaling.
- **G-vector/Band**: For large setups (single Gamma point), rely on G-vector and Band parallelization.

### Optimization:
- **Process Binding**: Disable hyperthreading for improved floating-point performance.
- **Memory**: Ensure sufficient memory per core; excessive swapping kills performance.
- **Start-up**: Use a converged wavefunction from a smaller cutoff or similar system to speed up SCF.

## Community and Support
- **Primary Hub**: [Psi-k Portal](https://psi-k.net/codes/petot).
- **Development**: Historically centered at LBNL (Lin-Wang Wang group).
- **Support Channels**: Direct contact with developers or academic collaboration; less community forum traffic than VASP.

## Limitations & Known Constraints
- **Availability**: Standard distribution is less centralized than codes like VASP or Quantum ESPRESSO; often obtained via direct academic contact or specific archives.
- **Documentation**: Less comprehensive or modern online documentation compared to major community codes.
- **User Interface**: Typically relies on traditional text-based input files without a modern GUI.

## Comparison with Other Codes
- **vs VASP/QE**: PEtot focuses heavily on large-scale parallel performance specifically for plane-wave calculations, though it lacks the vast feature set (e.g., advanced functionals, extensive post-processing) of VASP or QE.
- **vs BigDFT**: PEtot uses standard plane waves, whereas BigDFT uses wavelets. PEtot is a traditional plane-wave code optimized for size.

## Verification & Sources
**Primary sources**:
1. Psi-k Code Database: https://psi-k.net/codes/petot
2. L.-W. Wang Group Publications (LBNL)
3. "PEtot: A plane-wave pseudopotential density functional theory program for large systems" (Description in literature)

**Confidence**: VERIFIED - Code existence and primary features are well-documented in scientific literature and community databases.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: DFT/Plane-Wave
- Key Feature: Massive Parallelism
