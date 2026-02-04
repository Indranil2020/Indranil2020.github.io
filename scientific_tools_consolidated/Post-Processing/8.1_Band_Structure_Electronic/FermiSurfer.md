# FermiSurfer

## Official Resources
- Homepage: https://fermisurfer.osdn.jp/
- Documentation: https://fermisurfer.osdn.jp/en/_build/html/index.html
- Source Repository: https://github.com/mitsuaki1987/fermisurfer
- License: MIT License

## Overview
FermiSurfer is a visualization tool for Fermi surfaces. It can display Fermi surfaces with color plots of physical quantities (such as Fermi velocity, superconducting gap, spin character) on the surface. It is designed to be lightweight, fast, and easy to use, supporting input from major first-principles codes.

**Scientific domain**: Electronic structure visualization, Fermi surfaces  
**Target user community**: Condensed matter physicists, DFT users

## Capabilities (CRITICAL)
- **Fermi Surface Visualization**: 3D interactive plotting
- **Property Mapping**: Color mapping of scalar quantities on the Fermi surface
- **Stereoscopic View**: Anaglyph 3D support
- **Nodal Lines**: Highlighting nodal lines and points
- **Input Support**: Native formats for Quantum ESPRESSO, interacting with VASP via tools
- **Cross-Platform**: Linux, macOS, Windows

**Sources**: FermiSurfer documentation, Comp. Phys. Comm. 239, 272 (2019)

## Inputs & Outputs
- **Input formats**: .frm (FermiSurfer format), BXSF (XCrySDen format)
- **Output data types**: Screen capture, interactive view

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: `fs.x` post-processing tool
- **VASP**: Can be converted via `vaspkit` or `ifermi`
- **AkaiKKR**: Supported
- **OpenMX**: Supported

## Workflow and Usage
1. Perform SCF/NSCF calculation.
2. Generate Fermi surface data (e.g., using `fs.x` in QE).
3. Run FermiSurfer: `fermisurfer fermi_surface.frm`
4. Interact with the GUI to rotate, slice, and color-map.

## Performance Characteristics
- Lightweight and fast (OpenGL based)
- Handles complex Fermi surface topologies

## Application Areas
- Superconductivity (gap symmetry visualization)
- Magnetism (spin projection)
- Transport properties (Fermi velocity mapping)
- Topological materials (Weyl points, nodal lines)

## Community and Support
- Open-source (MIT)
- Developed by Mitsuaki Kawamura (ISSP, U. Tokyo)
- Active on GitHub

## Verification & Sources
**Primary sources**:
1. Homepage: https://fermisurfer.osdn.jp/
2. GitHub: https://github.com/mitsuaki1987/fermisurfer
3. Publication: M. Kawamura, Comp. Phys. Comm. 239, 272 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Kawamura)
- Applications: Fermi surface visualization, property mapping, QE integration
