# Artaios

## Official Resources
- **Homepage**: https://www.chemie.uni-hamburg.de/institutes/ac/research-groups/herrmann/software/artaios.html
- **License**: Free for Academic Use (Request required)

## Overview
**Artaios** is a specialized post-processing code designed for the analysis of **electron transport in molecular junctions**. It bridges the gap between quantum chemical calculations (from Gaussian, ORCA, or similar) and transport phenomenology. By utilizing the Non-Equilibrium Green's Function (NEGF) formalism—typically in the wide-band limit or with rigorous self-energies—Artaios computes transmission probabilities, current-voltage characteristics, and local current pathways to elucidate the mechanisms of charge flow through individual molecules.

**Scientific domain**: Molecular Electronics, Spintronics, Chemical Physics
**Target user community**: Chemists and physicists designing single-molecule devices

## Theoretical Methods
- **Post-Processing NEGF**: Takes the Fock/Kohn-Sham, Overlap, and self-energy matrices to construct the Green's function.
- **Landauer-Büttiker**: Calculates coherent transport observables.
- **Local Currents**: Decomposes the total current into bond-currents to visualize pathways (e.g., through $\pi$-systems vs $\sigma$-bonds).
- **Interference Analysis**: Tools to identify constructive/destructive quantum interference (QI) in molecular wires.

## Capabilities
- **Observables**:
  - Transmission function $T(E)$.
  - Current-Voltage (I-V) curves.
  - Seebeck coefficient (Thermoelectric power).
  - Spin polarization (for spintronics).
- **Analysis**:
  - COOP/COHP analysis adapted for transport (Bond currents).
  - Molecular orbital projection.

## Key Strengths
- **Chemical Intuition**: Provides tools (like local currents) that map transport physics back to chemical concepts (bonds, orbitals).
- **Flexibility**: Interfaced with popular quantum chemistry codes (Gaussian, ORCA, Q-Chem, ADF).
- **Efficiency**: As a post-processing tool, it avoids the heavy computational cost of running a full self-consistent NEGF-DFT loop for every geometry, making it ideal for screening large sets of conformers.

## Inputs & Outputs
- **Inputs**:
  - Output files from Quantum Chemistry codes (formatted checkpoint files).
  - `artaios.inp`: Control file.
- **Outputs**:
  - `transmission.dat`: T(E) data.
  - `current_map.xyz`: Vector field of local currents for visualization.

## Interfaces & Ecosystem
- **Upstream**: Gaussian, ORCA, Q-Chem, NWChem, ADF.
- **Visualization**: VMD or similar tools for viewing current maps.

## Performance Characteristics
- **Speed**: Very fast ($O(N^3)$ of the basis size, but basis is usually small for single molecules).
- **Scalability**: Embarrassingly parallel over energy points or different molecular conformers.

## Limitations & Known Constraints
- **Self-Consistency**: Does not re-optimize the density under bias (non-self-consistent); accurate only for low bias or linear response.
- **Electrostatics**: Often relies on effective screening models rather than solving Poisson equation.

## Comparison with Other Codes
- **vs. GOLLUM**: Similar post-processing philosophy; Artaios has strong features for "local current" visualization and chemical bond analysis.
- **vs. Smeagol/TranSIESTA**: These are fully self-consistent NEGF-DFT codes; Artaios is a lighter post-processing layer.

## Application Areas
- **Molecular Switches**: Designing molecules that change conductance upon photo-switching.
- **Quantum Interference**: Engineering molecules with Fano resonances.
- **Spintronics**: Spin-filtering in chiral molecules (CISS effect studies).

## Community and Support
- **Development**: Herrmann Group (University of Hamburg).
- **Source**: Distributed upon request from the group website.

## Verification & Sources
- **Website**: [https://www.chemie.uni-hamburg.de/...](https://www.chemie.uni-hamburg.de/institutes/ac/research-groups/herrmann/software/artaios.html)
- **Primary Publication**: C. Herrmann et al. (Check website).
- **Verification status**: ✅ VERIFIED
  - Established research code in molecular electronics.
