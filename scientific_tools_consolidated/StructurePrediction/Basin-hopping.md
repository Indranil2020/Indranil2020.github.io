# Basin hopping

## Official Resources
- Homepage: Method - Implemented in ASE, Wales Group codes, etc.
- Documentation: https://wiki.fysik.dtu.dk/ase/ase/optimize.html#basinhopping (ASE implementation)
- Source Repository: https://gitlab.com/ase/ase
- License: Varies (ASE is LGPL)

## Overview
Basin hopping is a global optimization algorithm used to find the global minimum of a potential energy surface. It transforms the energy landscape into a set of basins of attraction and samples them using Monte Carlo moves followed by local minimization. It is highly effective for cluster structure prediction and finding stable configurations of molecules and defects.

**Scientific domain**: Global optimization, structure prediction, clusters  
**Target user community**: Computational chemists, materials scientists

## Theoretical Methods
- Basin-Hopping (BH) / Monte Carlo Minimization (MCM)
- Canonical Monte Carlo sampling
- Local minimization (L-BFGS, CG, etc.)
- Metropolis criterion
- Temperature-dependent sampling

## Capabilities (CRITICAL)
- Global optimization of atomic structures
- Prediction of cluster geometries (Lennard-Jones, metallic clusters)
- Surface reconstruction search
- Adsorbate configuration search
- Implemented in: ASE, GMIN (Wales group), LAMMPS (via Python/fix), etc.

**Sources**: D.J. Wales and J.P.K. Doye, J. Phys. Chem. A 101, 5111 (1997)

## Inputs & Outputs
- **Input formats**: Initial structure, potential/calculator, temperature
- **Output data types**: Lowest energy structure, trajectory of minima

## Interfaces & Ecosystem
- **ASE**: Standard Python implementation (`ase.optimize.BasinHopping`)
- **Potentials**: Works with any calculator (DFT, classical)

## Performance Characteristics
- Stochastic, requires many minimization steps
- Embarrassingly parallel (multiple independent runs)
- Efficient for cluster physics

## Application Areas
- Nano-clusters
- Protein folding (simplified models)
- Defect searching
- Molecular docking

## Community and Support
- Broad usage in chemical physics
- ASE community
- Wales group (Cambridge)

## Verification & Sources
**Primary sources**:
1. Publication: D.J. Wales and J.P.K. Doye, J. Phys. Chem. A 101, 5111 (1997)
2. ASE Documentation: https://wiki.fysik.dtu.dk/ase/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Method: STANDARD
- Documentation: AVAILABLE (ASE)
- Applications: Global optimization, structure prediction, clusters
