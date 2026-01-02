# WannierTools

## Official Resources
- Homepage: http://www.wanniertools.com/
- Documentation: http://www.wanniertools.org/tutorials.html
- Source Repository: https://github.com/quanshengwu/wannier_tools
- License: GNU General Public License v3.0

## Overview
WannierTools is an open-source software package for investigating topological properties of materials using tight-binding models from Wannier functions. It calculates topological invariants, surface states, and various spectroscopic properties for topological materials.

**Scientific domain**: Topological materials, surface states, Berry phase physics  
**Target user community**: Researchers studying topological insulators, Weyl semimetals, and topological properties

## Theoretical Methods
- Tight-binding Hamiltonian from Wannier functions
- Berry curvature and Berry phase calculations
- Wilson loop methods
- Surface Green's function techniques
- k·p perturbation theory
- Symmetry analysis

## Capabilities (CRITICAL)
- Topological invariants (Z2, Chern number, mirror Chern number)
- Surface state calculations (surface spectrum, Fermi arcs)
- Berry curvature and anomalous Hall conductivity
- Orbital Hall conductivity
- Spin Hall conductivity
- Chiral anomaly calculations
- Landau level spectrum
- Weyl/Dirac point identification
- Nodal line calculations
- 3D Fermi surface plotting
- Energy band structures along arbitrary paths
- Joint density of states
- Optical conductivity
- Quantum metric tensor

**Sources**: Official WannierTools documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - wt.in (main input file)
  - wannier90_hr.dat (tight-binding Hamiltonian from Wannier90)
  - wannier90_centres.xyz (Wannier center positions)
  
- **Output data types**:
  - Surface state spectra
  - Topological invariant values
  - Berry curvature data
  - Fermi surface files
  - Band structure data
  - Conductivity tensors

## Interfaces & Ecosystem
- **Wannier90 integration**:
  - Direct reading of Wannier90 outputs
  - Tight-binding Hamiltonians from Wannier90
  
- **DFT interfaces**:
  - Via Wannier90 (VASP, Quantum ESPRESSO, WIEN2k, etc.)
  
- **Visualization**:
  - Gnuplot scripts generated
  - ParaView compatible outputs
  - Python plotting scripts

## Limitations & Known Constraints
- **Requires Wannier functions**: Needs converged Wannier90 calculation first
- **Tight-binding approximation**: Limited to TB model accuracy
- **Computational cost**: Surface state calculations can be expensive
- **Learning curve**: Requires understanding of topological concepts
- **Documentation**: Good but assumes topology knowledge
- **Platform**: Linux/Unix; Fortran compilation required

## Verification & Sources
**Primary sources**:
1. Official website: http://www.wanniertools.com/
2. Documentation: http://www.wanniertools.org/tutorials.html
3. GitHub repository: https://github.com/quanshengwu/wannier_tools
4. Q. Wu et al., Comput. Phys. Commun. 224, 405 (2018) - WannierTools paper

**Secondary sources**:
1. WannierTools tutorials and examples
2. Published topological materials studies
3. Wannier90 integration documentation
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (GitHub issues, email)
- Academic citations: >200
- Active development: Regular updates
