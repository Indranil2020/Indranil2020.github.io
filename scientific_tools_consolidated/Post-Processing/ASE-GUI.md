# ASE-GUI

## Official Resources
- Homepage: https://wiki.fysik.dtu.dk/ase/ase/gui/gui.html
- Documentation: https://wiki.fysik.dtu.dk/ase/ase/gui/gui.html
- Source Repository: https://gitlab.com/ase/ase
- License: GNU Lesser General Public License v2.1

## Overview
ASE-GUI is the graphical user interface for the Atomic Simulation Environment (ASE). It is a lightweight, Python-based viewer that allows users to visualize, manipulate, and analyze atomic structures. It is integrated directly into ASE and can handle any file format supported by ASE (xyz, cif, POSCAR, trajectory files, etc.). It is particularly useful for quick inspection of structures, setting up constraints, and visualizing MD trajectories or relaxation paths.

**Scientific domain**: Molecular visualization, structure manipulation, trajectory viewing  
**Target user community**: ASE users, computational materials scientists

## Theoretical Methods
- Visualization of atomic coordinates and bonds
- Crystal lattice visualization (unit cells)
- Trajectory animation
- Graph plotting (Energy vs time, Forces)
- Reciprocal space visualization (Brillouin zone)

## Capabilities (CRITICAL)
- **Visualization**: Atoms, unit cells, bonds, velocities
- **Structure Manipulation**: Move atoms, rotate/translate, add vacuum, repeat unit cell
- **Analysis**: Measure distances, angles, dihedrals
- **Trajectory Viewing**: Playback of MD or relaxation trajectories (.traj, .bundle)
- **Plotting**: Built-in plotting of energy, forces, and other properties from trajectory files
- **Surfaces**: Create surfaces and adsorbates (via ASE build tools)
- **Movie Export**: Save animations as GIF/MP4

**Sources**: ASE documentation

## Inputs & Outputs
- **Input formats**: Any format supported by `ase.io` (VASP, Gaussian, XYZ, CIF, PDB, etc.)
- **Output data types**: Images (PNG, EPS), Structure files, Movies

## Interfaces & Ecosystem
- **ASE**: Native part of ASE
- **Matplotlib**: Used for plotting graphs
- **Python**: Accessible via `ase gui` command or `view()` function in scripts

## Workflow and Usage
1. Open structure: `ase gui structure.xyz`
2. Open trajectory: `ase gui optimization.traj`
3. Manipulate: Select atoms, use tools to rotate or modify.
4. Graph: `Graph` -> `Energy` to see relaxation progress.
5. Save: `File` -> `Save` to export modified structure.

## Performance Characteristics
- Lightweight and fast for small to medium systems
- Not designed for rendering millions of atoms (use OVITO/VMD for that)
- Quick startup time

## Application Areas
- Checking convergence of geometry optimizations
- Validating initial structures
- Visualizing NEB paths
- Quick measurements of bond lengths

## Community and Support
- Part of ASE ecosystem
- Active mailing list and GitLab repository
- Developed by ASE community (DTU and contributors)

## Verification & Sources
**Primary sources**:
1. Documentation: https://wiki.fysik.dtu.dk/ase/ase/gui/gui.html
2. GitLab: https://gitlab.com/ase/ase

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitLab)
- Development: ACTIVE
- Applications: Visualization, trajectory analysis, structure editing
