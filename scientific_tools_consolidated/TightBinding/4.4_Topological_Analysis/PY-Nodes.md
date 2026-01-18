# PY-Nodes

## Official Resources
- **Homepage**: https://sourceforge.net/projects/py-nodes/
- **Repository**: https://sourceforge.net/projects/py-nodes/
- **License**: GPL-3.0

## Overview
**PY-Nodes** is a Python-based computational tool designed to automatically search for and classify **band degeneracy points** (nodes) in the Brillouin zone of topological semimetals. Specifically tailored for the all-electron DFT code **WIEN2k**, it uses a simplex optimization algorithm (Nelder-Mead) to minimize the energy gap function $\Delta E(\mathbf{k}) = |E_{n+1}(\mathbf{k}) - E_n(\mathbf{k})|$. This allows it to locate Weyl points, Dirac points, and nodal lines with high precision without computationally expensive dense grid scans.

**Scientific domain**: Topological Semimetals, Optimization
**Target user community**: WIEN2k users searching for Weyl points

## Theoretical Methods
- **Optimization**: Nelder-Mead method to find local minima of the gap function.
- **Topological Nodes**:
  - **Weyl Points**: Two-band crossings in 3D.
  - **Dirac Points**: Four-band crossings (typically protected by symmetry).
  - **Nodal Lines**: Continuous loops of degeneracy.
- **DFT backend**: WIEN2k (FLAPW method).

## Capabilities
- **Search**:
  - Automatic detection of gap closing points.
  - Can trace nodal lines by following the degeneracy valley.
- **Classification**:
  - Distinguishes between point nodes and lines based on the Hessian of the gap.
- **precision**: Finds coordinates to machine precision, unrestricted by a pre-defined k-mesh.

## Key Strengths
- **All-Electron Accuracy**: By using WIEN2k, it is suitable for f-electron systems and heavy metals where pseudopotential errors might shift node positions.
- **Efficiency**: Much faster than grid-based methods ($O(N_{iter})$ vs $O(N_k^3)$), essential for searching the full 3D BZ.
- **Automation**: Can be scripted to scan multiple band pairs.

## Inputs & Outputs
- **Inputs**:
  - WIEN2k `case.energy` files.
  - Search configuration (start points, bands).
- **Outputs**:
  - List of node coordinates ($k_x, k_y, k_z$) and residual gaps.

## Interfaces & Ecosystem
- **Upstream**: WIEN2k.
- **Dependencies**: NumPy.

## Performance Characteristics
- **Speed**: The search algorithm is very fast; the bottleneck is typically reading the initial DFT data or running the DFT steps if on-the-fly calculation is needed (optional workflow).
- **Reliability**: Requires reasonable starting guesses to avoid getting stuck in local (non-zero) minima.

## Comparison with Other Codes
- **vs. WannierTools**: WannierTools finds nodes using a tight-binding model. PY-Nodes works directly with the DFT eigenvalues, avoiding Wannierization errors but limited by DFT cost if dynamic recalculation is needed.
- **vs. IrRep**: IrRep finds nodes at *high symmetry points* via representations. PY-Nodes searches for accidental crossings at *generic* k-points (Weyl points).

## Application Areas
- **Weyl Semimetals**: TaAs, NbP.
- **Nodal Line Semimetals**: ZrSiS family.

## Community and Support
- **Development**: V. Pandey and S.K. Pandey.
- **Source**: SourceForge.

## Verification & Sources
- **Primary Publication**: V. Pandey et al., Comput. Phys. Comm. (2023).
- **Verification status**: ✅ VERIFIED
  - Published research code.
