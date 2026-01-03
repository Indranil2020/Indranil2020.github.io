# pymatgen (Python Materials Genomics)

## Official Resources
- Homepage: https://pymatgen.org/
- Documentation: https://pymatgen.org/
- Source Repository: https://github.com/materialsproject/pymatgen
- License: MIT License

## Overview
Pymatgen (Python Materials Genomics) is a robust, open-source Python library for materials analysis. It provides a core set of objects to represent materials (e.g., structures, molecules) and a comprehensive suite of tools to generate input files for and parse output files from various electronic structure codes. It is the core analysis code powering the **Materials Project**.

**Scientific domain**: Materials science, high-throughput analysis, crystallography
**Target user community**: Computational materials scientists, data scientists

## Capabilities
- **Structure Objects**: Robust classes for defining and manipulating crystal structures (`Structure`) and molecules (`Molecule`).
- **Input/Output**: Generation of input sets and parsing of outputs for major codes:
  - VASP, Quantum ESPRESSO, CP2K, ABINIT, CASTEP, LAMMPS, QChem, NWChem, etc.
- **Analysis**:
  - Phase diagrams (PD construction, stability analysis)
  - Electronic structure (Band structure, DOS plotting)
  - Diffusion analysis (NEB path finding, diffusion coefficient)
  - Diffraction patterns (XRD, neutron)
  - Pourbaix diagrams
- **Symmetry**: Integration with **spglib** for symmetry analysis.
- **Database**: Tools to interact with the Materials Project API (`MPRester`).

## Interfaces & Ecosystem
- **Integration**: Works seamlessly with **FireWorks** (workflows), **Custodian** (error handling), and **Atomate** (pre-built workflows).
- **Dependencies**: NumPy, SciPy, Matplotlib, Spglib, NetworkX, Pandas.

## Workflow and Usage
Pymatgen is typically used as a library within Python scripts or Jupyter notebooks.

```python
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

# Load structure from file
structure = Structure.from_file("POSCAR")

# Analyze symmetry
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
sga = SpacegroupAnalyzer(structure)
print(sga.get_space_group_symbol())

# Generate input files
poscar = Poscar(structure)
poscar.write_file("POSCAR_new")
```

## Application Areas
- High-throughput materials screening
- Phase stability analysis
- Electronic structure analysis
- Data generation for machine learning models

## Verification & Sources
**Primary sources**:
1. Homepage: https://pymatgen.org/
2. GitHub: https://github.com/materialsproject/pymatgen

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: High-throughput analysis, Materials Project core
