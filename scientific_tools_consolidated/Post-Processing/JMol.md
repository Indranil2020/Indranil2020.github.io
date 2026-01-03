# Jmol

## Official Resources
- Homepage: http://jmol.sourceforge.net/
- Documentation: http://jmol.sourceforge.net/docs/
- Source Repository: https://sourceforge.net/projects/jmol/
- License: GNU Lesser General Public License v2.1

## Overview
Jmol is a free, open-source Java viewer for chemical structures in 3D. It is designed to be cross-platform, running on Windows, Mac, Linux, and Unix systems. It supports a huge range of file formats and is widely used for teaching, research, and web-based visualization (via JSmol, the HTML5/JavaScript version). It can visualize molecules, crystals, materials, and biomolecules, along with orbitals, surfaces, and vibrations.

**Scientific domain**: Molecular visualization, crystallography, chemical education  
**Target user community**: Chemists, biochemists, materials scientists, educators

## Theoretical Methods
- 3D Rendering (Java 3D / WebGL for JSmol)
- Isosurface generation (Marching Cubes)
- Symmetry operations (Space groups)
- Vibrational animation
- Molecular orbital visualization (Gaussian/Slater)

## Capabilities (CRITICAL)
- **Visualization**: Atoms, bonds, polyhedra, surfaces, orbitals, vibrations
- **Web Integration**: JSmol allows 3D visualization in browsers without plugins
- **File Support**: Reads >60 formats (CIF, PDB, XYZ, Gaussian, VASP, QE, etc.)
- **Scripting**: Powerful scripting language for animation and interaction
- **Symmetry**: Automatic detection and display of symmetry elements
- **Measurement**: Distances, angles, torsion angles

**Sources**: Jmol website

## Inputs & Outputs
- **Input formats**: CIF, PDB, XYZ, CML, MOL, Gaussian logs, VASP POSCAR/OUTCAR, etc.
- **Output data types**: Images (JPG, PNG), POV-Ray scenes, VRML, script files

## Interfaces & Ecosystem
- **JSmol**: JavaScript version for web embedding
- **Java**: Native Java application
- **Browser**: Runs in any modern web browser (JSmol)

## Workflow and Usage
1. Open file: `File` -> `Open` (or drag and drop).
2. Manipulate view: Rotate, zoom, select atoms.
3. Scripting: Open console, type commands (e.g., `color atoms cpk`, `isosurface cutoff 0.02 mo 10`).
4. Export: `File` -> `Export` -> `Image`.

## Performance Characteristics
- Lightweight and fast for small/medium systems
- Java version faster than JSmol for large structures
- Highly portable

## Application Areas
- Chemical education (interactive tutorials)
- Database visualization (PDB, Crystallography Open Database)
- Quick inspection of calculation results
- Website integration

## Community and Support
- Open-source (LGPL)
- Long-standing project (>20 years)
- Active mailing list and wiki
- Developed by Bob Hanson and others

## Verification & Sources
**Primary sources**:
1. Homepage: http://jmol.sourceforge.net/
2. Wiki: http://wiki.jmol.org/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (SourceForge)
- Development: ACTIVE
- Applications: Molecular visualization, JSmol, education
