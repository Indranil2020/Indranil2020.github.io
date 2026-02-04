# PyMOL

## Official Resources
- Homepage: https://pymol.org/
- Documentation: https://pymolwiki.org/index.php/Main_Page
- Source Repository: https://github.com/schrodinger/pymol-open-source
- License: Open-source (Open source build) / Commercial (Incentive PyMOL)

## Overview
PyMOL is a molecular visualization system on an open-source foundation, maintained by Schrödinger. It is renowned for its high-quality rendering of biomolecules (proteins, nucleic acids) and its powerful Python scripting capabilities. PyMOL excels at creating publication-quality images and movies, and performing structural analysis like superposition and measurement.

**Scientific domain**: Molecular visualization, structural biology, rendering  
**Target user community**: Structural biologists, biochemists, computational chemists

## Theoretical Methods
- Ray tracing (built-in ray tracer)
- Molecular surface generation (solvent accessible surface)
- Electrostatic potential mapping (APBS plugin)
- Structural alignment (superposition algorithms)
- Molecular editing (mutagenesis, torsion adjustment)

## Capabilities (CRITICAL)
- **High-Quality Rendering**: Ray-tracing for publication figures (shadows, occlusion)
- **Representation**: Cartoon, ribbons, sticks, spheres, surfaces, mesh, dots
- **Analysis**: RMSD alignment, hydrogen bond detection, distance measurement
- **Scripting**: Full Python API (`pymol` module) and command language (`pml` scripts)
- **Session Saving**: Save state as `.pse` files
- **Plugins**: Extensive ecosystem (e.g., APBS for electrostatics, Autodock for docking)

**Sources**: PyMOL Wiki, Schrödinger website

## Inputs & Outputs
- **Input formats**: PDB, mmCIF, SDF, MOL2, XYZ, CCP4 map, XPLOR map, cube, DX
- **Output data types**: PNG/Ray-traced images, MPEG movies, PDB files, session files (.pse)

## Interfaces & Ecosystem
- **Python**: Deep integration (PyMOL IS a Python interpreter extended with C libraries)
- **Jupyter**: Can be run within Jupyter notebooks (via ipymol or similar)
- **APBS**: Integration for electrostatics
- **VMD**: Complementary tool (VMD often for dynamics, PyMOL for static high-quality images)

## Workflow and Usage
1. Load structure: `load protein.pdb`
2. Style: `hide everything`, `show cartoon`, `show sticks, resn LIG`
3. Color: `color red, chain A`, `color blue, chain B`
4. Ray trace: `ray 2000, 2000`
5. Save image: `png image.png`

## Performance Characteristics
- Highly optimized rendering engine
- Can handle large protein complexes
- Open-source build may lack some performance features of the commercial Incentive PyMOL

## Application Areas
- Structural biology (crystallography, Cryo-EM)
- Drug design (docking visualization)
- Molecular dynamics snapshots
- Educational illustrations

## Community and Support
- PyMOLWiki (community maintained)
- PyMOL-users mailing list
- Open source and commercial versions

## Verification & Sources
**Primary sources**:
1. Homepage: https://pymol.org/
2. Wiki: https://pymolwiki.org/
3. GitHub: https://github.com/schrodinger/pymol-open-source

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE (Wiki)
- Source: OPEN (GitHub)
- Development: ACTIVE (Schrödinger)
- Applications: High-quality rendering, structural biology, Python scripting
