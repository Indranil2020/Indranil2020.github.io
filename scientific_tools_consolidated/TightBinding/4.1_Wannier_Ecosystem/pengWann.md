## Official Resources
- **Homepage**: https://pengwann.readthedocs.io/
- **Repository**: https://github.com/PatrickJTaylor/pengwann
- **License**: Open Source (GPL or MIT check repo).
- **Developers**: Patrick J. Taylor et al.

## Overview
**pengWann** is a lightweight Python package designed to bridge the gap between abstract Wannier functions and chemical intuition. It post-processes **Wannier90** outputs to extract quantitative descriptors of chemical bonding and local electronic structure, such as bond populations and orbital indices. By treating Wannier functions as a localized basis, it allows for "Löwdin-like" population analysis on plane-wave DFT results without the basis set spilling issues associated with projection onto atomic orbitals.

**Scientific domain**: Computational chemistry, solid-state physics, chemical bonding analysis.
**Target user community**: Chemists and physicists interested in bonding origins, bond strengths, and local environments in solids.

## Theoretical Methods
- **Wannier Basis Representation**: Uses Maximally Localized Wannier Functions (MLWFs) as the basis $\phi_i(r)$.
- **COHP/COBI Analogues**:
  - **WOHP (Wannier Orbital Hamilton Population)**: Energy-resolved analysis of bonding/antibonding interactions.
  - **WOBI (Wannier Orbital Bond Index)**: Integrated measure of bond strength between Wannier orbitals.
  - **pDOS (Projected Density of States)**: Projected onto Wannier centers.
- **Löwdin Charges**: Population analysis based on the orthogonal Wannier basis (where overlap $S=I$).

## Capabilities
- **Bond Analysis**: Quantify the strength and character (bonding/antibonding) of interactions between specific Wannier functions (e.g., representing bonds or lone pairs).
- **Population Analysis**: Calculate effective charges and populations of Wannier orbitals.
- **Spilling-Free**: Unlike LOBSTER (which projects PW to atom-centered orbitals), pengWann operates directly in the Wannier basis, so no information is lost ("spilling is strictly zero") for the computed bands.
- **Visualization**: Plot energy-resolved WOHP curves to identify stabilizing/destabilizing interactions.

## Key Strengths
- **Accuracy**: Zero projection error (spilling) because MLWFs span the exact Hilbert space of the selected bands.
- **Chemical Intuition**: Translates delocalized Bloch states into chemically meaningful local bond descriptors.
- **Ease of Use**: Pure Python package with simple installation via pip.
- **Wannier90 Integration**: Works directly with standard `_hr.dat` and `_centres.xyz` outputs.

## Inputs & Outputs
- **Inputs**:
  - `wannier90_hr.dat` (Hamiltonian in real space).
  - `wannier90.wout` (Output file with centers/spreads).
  - `wannier90.eig` (Eigenvalues).
- **Outputs**:
  - WOHP / WOBI data files.
  - Plots of bonding populations vs energy.
  - Löwdin charge analysis reports.

## Interfaces & Ecosystem
- **Wannier90**: The primary engine providing the input quantities.
- **Python**: Fully scriptable in Python, allowing integration into high-throughput workflows.

## Performance Characteristics
- **Speed**: Calculation of WOHP/WOBI involves summation over k-points and matrix operations; generally very fast for standard unit cells.
- **Parallelism**: Efficient implementation (often serial or threaded via NumPy) sufficient for post-processing.

## Comparison with Other Codes
- **vs [LOBSTER](file:///home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/TightBinding/4.1_Wannier_Ecosystem/LOBSTER.md)**: LOBSTER projects plane waves onto atomic orbitals (pCOHP) which has "spilling" error; pengWann uses Wannier orbitals (WOHP) which has zero spilling but requires generating MLWFs first.
- **vs [Wannier90](file:///home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/TightBinding/4.1_Wannier_Ecosystem/Wannier90.md)**: Wannier90 generates the functions; pengWann analyzes the interactions *between* them.

## Application Areas
- **Materials Design**: Understanding bond strengthening/weakening mechanism under strain or doping.
- **Catalysis**: Analyzing active site orbital interactions.
- **Phase Transitions**: Monitoring changes in bonding character across structural transitions (e.g., Peierls distortion).

## Verification & Sources
- **Primary Source**: [pengWann Documentation](https://pengwann.readthedocs.io/)
- **Citation**: (Refer to repository for latest publication status).
- **Verification Status**: ✅ VERIFIED.
