# IrRep

## Official Resources
- **Homepage**: https://github.com/stepan-tsirkin/irrep
- **Documentation**: https://github.com/stepan-tsirkin/irrep/wiki
- **Repository**: https://github.com/stepan-tsirkin/irrep
- **License**: GPL-3.0

## Overview
**IrRep** is a Python code that calculates the **symmetry eigenvalues** (traces of representation matrices) of electronic Bloch states in crystals. By interfacing with standard DFT codes (VASP, Quantum ESPRESSO, Abinit, Wien2k), it allows researchers to identify the **irrep labels** of energy bands at high-symmetry points in the Brillouin Zone. This data is the primary input for the theory of **Topological Quantum Chemistry (TQC)**, enabling the diagnosis of topological phases via symmetry indicators.

**Scientific domain**: Topological Quantum Chemistry, Group Theory
**Target user community**: Researchers classifying topological materials

## Theoretical Methods
- **Group Theory**: 
  - Trace formulas for rotations/inversions acting on Bloch states.
  - Character tables for all 230 Space Groups (and magnetic groups).
- **Topological Quantum Chemistry**:
  - Symmetry Indicators ($\mathbb{Z}_2, \mathbb{Z}_4, \mathbb{Z}_{8}, \mathbb{Z}_{12}$).
  - Elementary Band Representations (EBRs) and compatibility relations.

## Capabilities
- **Calculation**:
  - Traces of Cuillin-Wigner Seitz operators.
  - Identification of Irreps (e.g., $\Gamma_6^+$, $X_5^-$).
- **Interfaces**:
  - **VASP** (`WAVECAR`).
  - **Quantum ESPRESSO** (`save` folder).
  - **Abinit** (`WFK`).
  - **Wien2k** (`case.energy`/`case.vector`).
- **Scope**:
  - Spinful (SOC) and spinless systems.
  - Non-magnetic and Magnetic Space Groups.

## Key Strengths
- **Universality**: Works with almost all major DFT codes, effectively standardizing the output for topological analysis.
- **Robustness**: Uses `spglib` for symmetry detection and handles non-symmorphic groups and double groups correctly, which is mathematically nontrivial.
- **Efficiency**: Only requires wavefunctions at a few high-symmetry k-points (TRIM points), making it orders of magnitude faster than calculating Wilson loops over the full BZ.

## Inputs & Outputs
- **Inputs**:
  - DFT wavefunction file.
  - Structure file (POSCAR).
  - K-points list.
- **Outputs**:
  - Tables of traces and irrep labels.
  - Files compatible with the Bilbao Crystallographic Server (BCS) tools.

## Interfaces & Ecosystem
- **BCS**: Output can be pasted into "CheckTop" or "Trace" tools on the Bilbao server.
- **WannierBerri**: Integrated for symmetry analysis in Wannier interpolation.

## Performance Characteristics
- **Speed**: Very fast (seconds).
- **scalability**: Limited only by the memory to load the DFT wavefunction of one k-point.

## Comparison with Other Codes
- **vs. Z2Pack**: Z2Pack calculates topology by "brute force" integration (Wilson loops). IrRep uses symmetry principles (TQC). IrRep is faster but applies only to phases protected by crystal symmetry (or inversion/TR).
- **vs. TopMat**: TopMat handles the classification logic for magnetic groups; IrRep is the engine that calculates the raw traces.

## Application Areas
- **Topological Semimetals**: Identifying band crossings protected by different irreps.
- **High-Throughput Screening**: Used in workflows to classify thousands of materials in databases (e.g., Materiae, Materials Project).

## Community and Support
- **Development**: Stepan Tsirkin (University of Zurich).
- **Source**: GitHub.

## Verification & Sources
- **Primary Publication**: I. Iraola et al., Comput. Phys. Comm. 272, 108226 (2022).
- **Verification status**: ✅ VERIFIED
  - Standard tool in the field.
