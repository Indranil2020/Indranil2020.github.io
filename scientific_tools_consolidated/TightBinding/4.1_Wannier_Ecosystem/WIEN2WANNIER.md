# WIEN2WANNIER

## Official Resources
- **Homepage**: https://wien2wannier.github.io/
- **Documentation**: https://github.com/wien2wannier/wien2wannier/wiki
- **Source Repository**: https://github.com/wien2wannier/wien2wannier
- **License**: GNU General Public License (GPL)

## Overview
**WIEN2WANNIER** is a specialized interface program that connects the high-precision, all-electron full-potential linearized augmented plane-wave (FP-LAPW) code **WIEN2k** with the maximally-localized Wannier function code **Wannier90**. It computes the necessary overlap matrices ($M_{mn}$) and projection matrices ($A_{mn}$) from WIEN2k's Bloch states, allowing for the construction of Wannier functions with all-electron accuracy. This tool is essential for researchers using methods that require a localized basis set (like DMFT or Berry Phase calculations) but demand the precision of LAPW.

**Scientific domain**: Condensed Matter Physics, Strongly Correlated Systems, Topological Materials
**Target user community**: Users of WIEN2k needing Wannier functions for post-processing

## Theoretical Methods
- **FP-LAPW Basis**: Utilizes the highly accurate all-electron wavefunctions from WIEN2k (Planewaves + Atomic spheres).
- **Projection**: Projects Bloch states onto a set of trial localized orbitals (s, p, d, f types) to generate initial guesses ($A_{mn}$).
- **Overlap Calculation**: Computes the overlaps between periodic parts of Bloch functions ($M_{mn}$) on a Monkhorst-Pack mesh.
- **Spin-Orbit Coupling**: Full support for non-collinear and spin-orbit coupled calculations via spinor projections.

## Capabilities
- **Wannier90 Interface**:
  - Generates `.mmn` (overlaps), `.amn` (projections), and `.eig` (eigenvalues) files.
  - Supports disentanglement of entangled bands.
- **All-Electron Accuracy**: Captures core and semi-core effects crucial for d- and f-electron systems.
- **Symmetry Handling**: Uses WIEN2k symmetry operations to reduce computational cost.
- **Real-Space Plotting**: Tools to visualize Wannier functions constructed from the LAPW basis (`wplot`).

## Key Strengths
- **Accuracy**: As an interface to WIEN2k, it provides MLWFs based on the "gold standard" of DFT methods (FP-LAPW).
- **f-electrons**: Particularly strong for lanthanides and actinides where pseudopotentials may struggle.
- **Integration**: Tightly integrated into the WIEN2k workflow (callable via `x w2w`).

## Inputs & Outputs
- **Inputs**:
  - WIEN2k structure (`case.struct`) and vector files (`case.vector`).
  - WIEN2WANNIER input file (`case.inw2w`).
  - Wannier90 input (`case.win`).
- **Outputs**:
  - Wannier90 required files: `case.mmn`, `case.amn`, `case.eig`.
  - Visualization data: `case.psZK`, `case.ploteig`.

## Interfaces & Ecosystem
- **WIEN2k**: Requires a working installation of WIEN2k.
- **Wannier90**: Generates inputs compatible with all modern versions of Wannier90.
- **dmft_proj**: Often used in conjunction with DMFT codes that interface with WIEN2k.

## Performance Characteristics
- **Computational Cost**: The interface step itself is relatively inexpensive compared to the SCF cycle; scaling depends on the number of k-points and bands.
- **Parallelism**: Supports k-point parallelization consistent with WIEN2k's MPI scheme.

## Limitations & Known Constraints
- **Complexity**: The LAPW basis is more complex than plane-waves, making the projection definition slightly more involved.
- **Dependencies**: Strictly tied to WIEN2k; cannot be used with other DFT codes.
- **Memory**: Overlap calculations for large systems/dense k-meshes can be memory-intensive.

## Comparison with Other Codes
- **vs. VASP2Wannier90**: VASP uses PAW potentials; WIEN2WANNIER uses all-electron LAPW. WIEN2WANNIER is preferred for systems where core states or orthogonality orthogonality is critical (e.g., NMR, hyperfine parameters, heavy fermions).
- **vs. SCAALD**: Another interface for all-electron codes, but WIEN2WANNIER is the official/standard one for WIEN2k.

## Application Areas
- **Strongly Correlated Materials**: DMFT studies of f-electron systems (Ce, Pu, U compounds).
- **Topological Insulators**: Accurate calculation of surface states and Berry curvature with SOC.
- **Fermi Surface Analysis**: High-precision Fermi surface interpolation via Wannier90.

## Community and Support
- **Development**: Maintained by the Institute of Solid State Physics, TU Wien (Elias Assmann, Peter Blaha et al.).
- **Mailing List**: Support provided via the active WIEN2k mailing list.
- **Updates**: Regularly updated with WIEN2k releases (e.g., v2.0 in WIEN2k 16.1).

## Verification & Sources
- **Official Website**: [https://wien2wannier.github.io/](https://wien2wannier.github.io/)
- **Primary Publication**: J. Kuneš et al., Comp. Phys. Commun. 181, 1888 (2010).
- **Verification status**: ✅ VERIFIED
  - Integral part of the WIEN2k suite.
  - Validated in numerous studies on f-electron dynamics.
