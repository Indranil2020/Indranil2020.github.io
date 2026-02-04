# Avogadro

## Official Resources
- Homepage: https://avogadro.cc/
- Documentation: https://avogadro.cc/docs/
- Source Repository: https://github.com/openchemistry/avogadro2
- License: GNU General Public License v2.0

## Overview
Avogadro is an advanced molecule editor and visualizer designed for cross-platform use in computational chemistry, molecular modeling, bioinformatics, materials science, and related areas. It offers flexible high quality rendering and a powerful plugin architecture. It is widely used for building molecular structures, generating input files for various quantum chemistry codes, and visualizing results.

**Scientific domain**: Molecular editing, visualization, input generation  
**Target user community**: Chemists, materials scientists, students

## Theoretical Methods
- Molecular Mechanics (MM) minimization (UFF, GAFF, MMFF94)
- Conformational searching
- Rendering (OpenGL)
- Isosurface generation (Orbitals, electrostatic potentials)
- Vibrational mode animation

## Capabilities (CRITICAL)
- **Structure Building**: Intuitive 3D builder for molecules, crystals, nanotubes, DNA/RNA
- **Input Generation**: Generates inputs for Gaussian, GAMESS, ORCA, Q-Chem, NWChem, VASP, LAMMPS, etc.
- **Visualization**: Orbitals, electron density surfaces, electrostatic potentials, vibrations
- **Analysis**: Bond lengths/angles, properties
- **Optimization**: Built-in geometry optimization using OpenBabel force fields
- **Extensibility**: Python scripting and plugin system

**Sources**: Avogadro website, J. Cheminform. 4, 17 (2012)

## Inputs & Outputs
- **Input formats**: CML, XYZ, PDB, Mol2, SDF, Gaussian/GAMESS output, Cube files
- **Output data types**: Structure files, images (PNG/JPG), POV-Ray scenes

## Interfaces & Ecosystem
- **OpenBabel**: Uses OpenBabel for file import/export and force fields
- **Quantum Codes**: Bridges GUI to HPC codes (input generation)
- **XtalOpt**: Hosts the XtalOpt structure prediction plugin

## Workflow and Usage
1. Draw molecule or import structure.
2. Optimize geometry: `Extensions` -> `Optimize Geometry`.
3. Setup calculation: `Extensions` -> `Gaussian` -> `Generate Input`.
4. Visualize results: Open output file, add `Surfaces` or `Vibrations`.

## Performance Characteristics
- Fast interactive rendering
- Efficient for small to medium molecules and unit cells
- Can struggle with massive MD trajectories compared to VMD

## Application Areas
- Preparation of DFT calculations
- Teaching chemistry
- Visualization of molecular orbitals
- Crystal structure building

## Community and Support
- Open-source (GPL)
- Part of OpenChemistry project
- Active development (Avogadro 2)
- Large user base

## Verification & Sources
**Primary sources**:
1. Homepage: https://avogadro.cc/
2. GitHub: https://github.com/openchemistry/avogadro2
3. Publication: M. D. Hanwell et al., J. Cheminform. 4, 17 (2012)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Kitware/OpenChemistry)
- Applications: Molecular editing, visualization, input generation
