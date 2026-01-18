# nested_wloop

## Official Resources
- **Homepage**: https://github.com/kuansenlin/nested_and_spin_resolved_Wilson_loop
- **Repository**: https://github.com/kuansenlin/nested_and_spin_resolved_Wilson_loop
- **License**: MIT License

## Overview
**nested_wloop** is a specialized Python toolkit designed to extend **PythTB** for the characterization of **Higher-Order Topological Insulators (HOTIs)** and **Fragile Topology**. It implements the numerical calculation of **Nested Wilson Loops**—a hierarchical Berry phase technique required to identify quadrupole and octupole insulators—as well as Spin-Resolved Wilson loops for systems with approximate time-reversal symmetry.

**Scientific domain**: Higher-Order Topology, Fragile Phases
**Target user community**: Theorists characterizing HOTIs and Twisted Bilayers

## Theoretical Methods
- **Nested Wilson Loop**: Diagonalization of the standard Wilson loop unitary $W_1$, followed by a second Wilson loop calculation along the remaining direction for the Wannier bands of $W_1$.
- **Spin-Resolved Wilson Loop**: Projecting the Wilson loop onto spin sectors ($P_{\uparrow} W P_{\uparrow}$) to define Spin Chern numbers in the absence of $S_z$ conservation.
- **Multipole Moments**: Calculation of bulk quadrupole moments $q_{xy}$ via nested geometric phases.

## Capabilities
- **HOTI Diagnosis**:
  - Identifies 2nd order topological insulators (corner states).
  - Calculates quadrupole invariants.
- **Fragile Topology**:
  - Detects "Fragile" bands that have trivial Chern number but non-trivial spin winding.
- **Workflow**:
  - Takes a `pythtb.tb_model` as input.
  - returns nested Wannier charge center spectra.

## Key Strengths
- **Cutting Edge**: One of the few public implementations of the *Benalcazar-Bernevig-Hughes (BBH)* nested loop procedure.
- **Integration**: Works directly with PythTB models, meaning users don't need to rewrite their Hamiltonians to use this advanced analysis.
- **Fragile Phases**: Essential for the modern study of "Wannier Obstructions" in twisted bilayer graphene and similar moiré systems.

## Inputs & Outputs
- **Inputs**: PythTB model object.
- **Outputs**: Plots of Nested Wannier Spectra ($\nu_y$ vs $k_x$).

## Interfaces & Ecosystem
- **Dependency**: PythTB (must be installed).
- **Language**: Python 3.

## Performance Characteristics
- **Complexity**: $O(N_k^2)$ due to the nested integration. Slower than standard Berry phase but manageable for typical tight-binding grids ($100 \times 100$).
- **Scalability**: Serial execution.

## Comparison with Other Codes
- **vs. Z2Pack**: Z2Pack does standard Wilson loops perfectly. nested_wloop adds the specific *nested* functionality for HOTIs which Z2Pack doesn't natively expose in a "one-shot" function.
- **vs. WannierTools**: WannierTools computes corner states via creating finite clusters (real space). nested_wloop predicts them from the *bulk* Bloch functions (momentum space).

## Application Areas
- **Quadrupole Insulators**: Checking the $\mathbb{Z}_2 \times \mathbb{Z}_2$ classification.
- **Twisted Bilayers**: Diagnosing fragile topology in flat bands.

## Community and Support
- **Development**: Kuan-Sen Lin.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/kuansenlin/nested_and_spin_resolved_Wilson_loop](https://github.com/kuansenlin/nested_and_spin_resolved_Wilson_loop)
- **Verification status**: ✅ VERIFIED
  - Specialized research extension.
