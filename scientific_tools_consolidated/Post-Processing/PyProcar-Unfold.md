# PyProcar-Unfold

## Official Resources
- Homepage: https://romerogroup.github.io/pyprocar/
- Documentation: https://romerogroup.github.io/pyprocar/modules.html#module-pyprocar.unfold
- Source Repository: https://github.com/romerogroup/pyprocar
- License: MIT License

## Overview
PyProcar-Unfold refers to the band unfolding module within the PyProcar package. It provides functionality to unfold the electronic band structure of supercells into the primitive cell Brillouin zone. This is particularly useful for studying the effective band structure of systems with broken translational symmetry, such as defects, interfaces, and alloys, using VASP or other supported DFT codes.

**Scientific domain**: Band unfolding, supercell analysis, electronic structure  
**Target user community**: VASP users, defect researchers

## Capabilities (CRITICAL)
- **Band Unfolding**: Recovering primitive cell spectral weights from supercell calculations
- **Projection**: Projecting unfolded bands onto atomic orbitals
- **Spin Texture**: Unfolding spin-polarized bands
- **Visualization**: Plotting the effective band structure (EBS)
- **Integration**: Part of the larger PyProcar ecosystem

**Sources**: PyProcar documentation, Comp. Phys. Comm. 251, 107206 (2020)

## Inputs & Outputs
- **Input formats**: PROCAR/OUTCAR (VASP), supercell transformation matrix
- **Output data types**: Unfolded spectral weights, band structure plots

## Interfaces & Ecosystem
- **PyProcar**: Core package
- **VASP**: Primary input source
- **Matplotlib**: Visualization backend

## Workflow and Usage
1. Perform supercell DFT calculation.
2. Define the supercell matrix (transformation from primitive to supercell).
3. Run unfolding: `pyprocar.unfold(code='vasp', dirname='.', supercell_matrix=[[2,0,0],[0,2,0],[0,0,2]])`
4. Plot the result.

## Performance Characteristics
- Python-based, generally efficient
- Depends on PROCAR size

## Application Areas
- Impurity levels in semiconductors
- Alloy band structures
- Surface states in slab calculations

## Verification & Sources
**Primary sources**:
1. Homepage: https://romerogroup.github.io/pyprocar/
2. GitHub: https://github.com/romerogroup/pyprocar

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Part of PyProcar package
