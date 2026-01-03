# SeeK-path

## Official Resources
- Homepage: https://seekpath.readthedocs.io/
- Documentation: https://seekpath.readthedocs.io/en/latest/
- Source Repository: https://github.com/giovannipizzi/seekpath
- License: MIT License

## Overview
SeeK-path is a python module and online tool for obtaining the standardized primitive cell and high-symmetry k-points path for band structure calculations. It automatically detects the space group of a crystal structure, converts it to a standard representation, and suggests a path through the Brillouin zone that captures all relevant features of the electronic bands.

**Scientific domain**: Band structure analysis, crystallography, symmetry analysis  
**Target user community**: DFT users, materials scientists

## Theoretical Methods
- Space group determination (via spglib)
- Standardization of crystal structures (conventions by HPkot, etc.)
- Brillouin zone analysis
- High-symmetry path generation
- Crystallographic conventions

## Capabilities (CRITICAL)
- Automatic detection of crystal symmetry
- Generation of standardized primitive cells
- Definition of high-symmetry k-path for band structures
- Web interface for quick visualization and path generation
- Python API for integration into workflows (e.g., AiiDA, ASE)
- Support for all 230 space groups

**Sources**: SeeK-path documentation, Comp. Mater. Sci. 128, 140 (2017)

## Inputs & Outputs
- **Input formats**: Crystal structure (POSCAR, CIF, Python objects)
- **Output data types**: Standardized structure, k-point coordinates, k-path labels, Brillouin zone geometry

## Interfaces & Ecosystem
- **AiiDA**: Deeply integrated into AiiDA workflows
- **Quantum ESPRESSO / VASP**: Generates KPOINTS inputs
- **Materials Cloud**: Powers the online SeeK-path tool
- **ASE / Pymatgen**: Compatible

## Workflow and Usage
1. Load structure: `structure = (cell, positions, numbers)`
2. Run SeeK-path: `res = seekpath.get_explicit_k_path(structure)`
3. Use output: `res['explicit_kpoints_rel']` for DFT input
4. Visualize: Use the generated k-path labels for plotting

## Performance Characteristics
- Fast symmetry analysis (spglib backend)
- Instantaneous for typical unit cells

## Application Areas
- High-throughput band structure calculations
- Database generation (standardized paths)
- Educational visualization of Brillouin zones

## Community and Support
- Open-source (MIT)
- Developed by Giovanni Pizzi (EPFL)
- Integrated into Materials Cloud

## Verification & Sources
**Primary sources**:
1. Homepage: https://seekpath.readthedocs.io/
2. GitHub: https://github.com/giovannipizzi/seekpath
3. Publication: Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Comp. Mater. Sci. 128, 140 (2017)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (EPFL/Materials Cloud)
- Applications: Standardized k-paths, symmetry, band structure workflows
