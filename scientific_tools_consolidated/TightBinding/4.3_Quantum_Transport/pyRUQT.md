# pyRUQT (Python Rowan University Quantum Transport)

## Official Resources
- **Repository**: https://github.com/HoyLab-Rowan/pyRUQT
- **License**: MIT License

## Overview
**pyRUQT** is a modular Python-based code for performing **Multiconfigurational Non-Equilibrium Green's Function (NEGF)** transport calculations. Unlike standard DFT-NEGF codes, pyRUQT integrates multiconfigurational electronic structure methods (like CASSCF/NEVPT2 or MC-PDFT) into the transport kernel. This allows for the accurate description of **strongly correlated** molecular junctions, capturing phenomena like the Coulomb blockade and Kondo effect that are often missed or poorly described by single-determinant DFT.

**Scientific domain**: Molecular Electronics, Strongly Correlated Systems, Quantum Transport
**Target user community**: Researchers in single-molecule transport and correlated electron physics

## Theoretical Methods
- **Multiconfigurational NEGF**: Combines NEGF with multiconfigurational wavefunctions.
- **MC-PDFT**: Multi-Configuration Pair-Density Functional Theory for efficient correlation.
- **Landauer-Büttiker**: Calculates transmission probabilities through the correlated region.
- **Coulomb Blockade**: Captures charge quantization and blockade steps.

## Capabilities
- **Transport Visualization**:
  - Transmission functions $T(E)$.
  - Current-Voltage (I-V) characteristics.
- **Correlation**:
  - Treats static and dynamic correlation in the scattering region.
  - Handles open-shell molecules and transition metal complexes.
- **Modularity**:
  - Driver-based design: Can interface with various electronic structure backends (e.g., PySCF, Molpro).

## Key Strengths
- **Accuracy**: Superior to DFT-NEGF for systems with near-degenerate states or localized spins.
- **Flexibility**: Python modularity allows easy testing of new functionals or hybrid schemes.
- **Accessibility**: Open-source and built on modern Python scientific stack.

## Inputs & Outputs
- **Inputs**:
  - Molecular geometry (XYZ).
  - Electronic structure parameters (from backend codes).
  - Transport setup (bias window, energy grid).
- **Outputs**:
  - `transmission.dat`: Energy-dependent transmission.
  - `current.dat`: Integrated current at different biases.

## Interfaces & Ecosystem
- **Upstream**:
  - **PySCF**: Primary open-source backend for generating Hamiltonian/Overlap matrices.
  - **Molpro**: Supported via file interface.
- **Downstream**:
  - **Matplotlib**: Visualization of transport curves.

## Performance Characteristics
- **Computational Cost**: dominated by the electronic structure calculation (CASSCF is expensive, scaling exponentially with active space size). The NEGF step is algebraic and relatively fast.
- **Scalability**: Parallelized over energy points; limited by the backend solver's scalability.

## Limitations & Known Constraints
- **Active Space size**: Limited to small active spaces ( ~16 orbitals) due to CASSCF complexity.
- **Self-Consistency**: Current implementation typically uses a "one-shot" or perturbative approach rather than fully self-consistent NEGF-SCF for the leads-molecule coupling.

## Comparison with Other Codes
- **vs. TranSIESTA / Smeagol**: These use DFT (Smeagol adds Hubbard U); pyRUQT uses explicit multiconfigurational wavefunctions, offering higher accuracy for strong correlation but at much higher cost.
- **vs. gDFTB**: Uses tight-binding/DFT; pyRUQT is for high-fidelity quantum chemistry transport.

## Application Areas
- **Single-Molecule Magnets**: Transport through magnetic molecules.
- **Radical Bridges**: Junctions with unpaired electrons.
- **Interference**: Destructive quantum interference in molecular wires.

## Community and Support
- **Development**: Developed by the Hoy Lab at Rowan University.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/HoyLab-Rowan/pyRUQT](https://github.com/HoyLab-Rowan/pyRUQT)
- **Primary Publication**: Garner et al., J. Chem. Phys. (Check repo for exact citation).
- **Verification status**: ✅ VERIFIED
  - Active academic project.
  - Methodology verified in literature.
