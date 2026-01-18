# TB2J

## Official Resources
- **Homepage**: https://tb2j.readthedocs.io/
- **Documentation**: https://tb2j.readthedocs.io/en/latest/
- **Source Repository**: https://github.com/mailhexu/TB2J
- **License**: BSD-2-Clause

## Overview
**TB2J** is an open-source Python package designed to calculate magnetic interaction parameters—specifically the isotropic Heisenberg exchange ($J$), anisotropic exchange, and Dzyaloshinskii-Moriya interaction (DMI)—from first-principles density functional theory (DFT) data. It utilizes the **Green's function method** combined with the **Magnetic Force Theorem**, allowing for the extraction of these parameters from a single unit-cell calculation without the need for computationally expensive supercell approaches. TB2J operates within a localized basis set framework, supporting both Maximally Localized Wannier Functions (MLWFs) and Linear Combination of Atomic Orbitals (LCAO).

**Scientific domain**: Magnetism, Spintronics, Topological Magnets, 2D Materials
**Target user community**: Researchers studying magnetic materials, magnon transport, and spin dynamics

## Theoretical Methods
- **Magnetic Force Theorem**: Treats small rotation of spins as a perturbation to the total energy.
- **Green's Function Approach**: Calculates inter-site magnetic interactions by integrating the Green's function, avoiding supercell total energy differences.
- **Localized Basis Sets**:
  - **Wannier Functions**: Via Wannier90 (for plane-wave codes like VASP, QE).
  - **LCAO**: Native support for numerical atomic orbitals (SIESTA, OpenMX, ABACUS).
- **Heisenberg Hamiltonian**: Maps DFT energetics to:
  $$H = -\sum_{i<j} J_{ij} \mathbf{S}_i \cdot \mathbf{S}_j + \sum_{i<j} \mathbf{D}_{ij} \cdot (\mathbf{S}_i \times \mathbf{S}_j) + \dots$$

## Capabilities
- **Interaction Parameters**:
  - Isotropic exchange ($J$).
  - Antisymmetric exchange / Dzyaloshinskii-Moriya Interaction (DMI, $\mathbf{D}$).
  - Anisotropic exchange tensor ($\Gamma$).
  - Biquadratic exchange (experimental).
- **Post-Processing**:
  - Magnon band structure calculation.
  - Spin-wave stiffness.
  - Curie temperature estimation (Mean Field Approx).
- **Data Generation**:
  - Input files for spin dynamics codes (Vampire, Multibinit, UppASD).
  - Visualization of magnetic interactions.

## Key Strengths
- **Efficiency**: Computes full distance-dependent interactions from a single unit-cell ground state calculation.
- **Versatility**: Works with almost all major DFT codes via Wannier90 or native interfaces.
- **Spin-Orbit Coupling**: Fully supports non-collinear calculations for extracting DMI and anisotropy.
- **Automation**: User-friendly Python interface simplifies the workflow from DFT to magnetic modeling.

## Inputs & Outputs
- **Inputs**:
  - Wannier90 outputs (`*_hr.dat`, `*.win`, `eig` files).
  - DFT Hamiltonian/Overlap matrices (for SIESTA/OpenMX/ABACUS).
  - Atomic structure and magnetic moments.
- **Outputs**:
  - Text files with $J_{ij}$ and $\mathbf{D}_{ij}$ lists.
  - `exchange.xml` (standard format).
  - Input files for Multibinit, Vampire, UppASD.
  - Plotting scripts for interactions.

## Interfaces & Ecosystem
- **Native Interfaces**:
  - **SIESTA**: Reads `.HS` and `.DE` files directly (`siesta2J.py`).
  - **OpenMX**: Reads `.scfout` and Hamiltonian files (`openmx2J.py`).
  - **ABACUS**: Uses LCAO output files (`abacus2J.py`).
- **Wannier90 Interface**:
  - Supports **VASP**, **Quantum ESPRESSO**, **Abinit**, **HK** and any code compatible with Wannier90 (`wannier90J.py`).
- **Downstream Integration**:
  - **Vampire**: Atomistic spin dynamics.
  - **Multibinit**: (Abinit project) Lattice dynamics and spin dynamics.
  - **UppASD**: Atomistic spin dynamics.

## Performance Characteristics
- **Speed**: The Green's function integration is lightweight; the dominant cost is the preceding DFT/Wannierization step.
- **Scaling**: Efficient for standard unit cells; cost scales with the number of basis functions (orbitals).
- **Parallelization**: Python multiprocessing supported for k-point integration.

## Limitations & Known Constraints
- **Mott Insulators**: Requires a valid DFT ground state; for strongly correlated systems (DFT+U), accurate U values are critical.
- **Metals**: RKKY interactions are captured, but convergence with k-points must be carefully checked.
- **Basis Quality**: Accuracy depends entirely on the quality of the Wannierization or the LCAO basis set.

## Comparison with Other Codes
- **vs. Total Energy Mapping**: Standard approach requires many large supercells to fit J values; TB2J uses a single cell and is much faster/less error-prone.
- **vs. Lichtenstein formula codes**: TB2J implements a variant of the Lichtenstein formula but generalizes to non-orthogonal/Wannier bases and DMI.
- **vs. SPR-KKR**: KKR methods also use Green's functions but are all-electron and muffin-tin based; TB2J works with standard pseudopotential codes.

## Application Areas
- **2D Magnets**: CrI3, Fe3GeTe2, etc. - studying distance-dependent exchange and thickness dependence.
- **Topological Magnets**: Skyrmion-hosting materials (DMI calculation).
- **Spintronics**: Antiferromagnets and complex magnetic textures.

## Community and Support
- **Source**: Hosted on GitHub.
- **Documentation**: ReadTheDocs.
- **Forum**: Issues and discussions on GitHub.
- **Tutorials**: Examples provided for all supported interfaces.

## Verification & Sources
- **Official Website**: [https://tb2j.readthedocs.io/](https://tb2j.readthedocs.io/)
- **Repository**: [https://github.com/mailhexu/TB2J](https://github.com/mailhexu/TB2J)
- **Primary Publication**: X. He et al., Computer Physics Communications 264, 107938 (2021).
- **Verification status**: ✅ VERIFIED
  - Active development (v0.7+).
  - Widely cited and used in magnetic materials research.
