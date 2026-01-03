# Pyiron

## Official Resources
- Homepage: https://pyiron.org/
- Documentation: https://pyiron.readthedocs.io/
- Source Repository: https://github.com/pyiron/pyiron
- License: BSD 3-Clause License

## Overview
Pyiron is an integrated development environment (IDE) for computational materials science. It provides a framework to manage the full lifecycle of simulations: from setting up structures and submitting jobs to analyzing data. It uses a project-based approach (HDF5 storage) and integrates with Jupyter notebooks to provide an interactive and reproducible workflow environment.

**Scientific domain**: Integrated Computational Materials Engineering (ICME), workflow management  
**Target user community**: Materials scientists, atomistic simulation users

## Capabilities (CRITICAL)
- **Unified Interface**: Wrappers for LAMMPS, VASP, SPHInX, DFTB+, and more.
- **Project Management**: organizing calculations in a hierarchical project structure.
- **HDF5 Storage**: Efficient storage of inputs and outputs in HDF5 format.
- **Interactive**: Designed to be used within Jupyter Notebooks.
- **Workflows**: Support for complex protocols (e.g., melting point, phase diagrams).
- **Structures**: Atomistic structure manipulation (based on ASE).

**Sources**: Pyiron website, Comp. Mater. Sci. 163, 24 (2019)

## Inputs & Outputs
- **Input formats**: Python objects
- **Output data types**: HDF5 files (`.h5`), generic output objects

## Interfaces & Ecosystem
- **Jupyter**: Primary user interface
- **ASE**: Used for structure handling
- **LAMMPS/VASP**: Core engines supported
- **pyiron_atomistics**: Module for atomistic simulation

## Workflow and Usage
1. Import: `from pyiron import Project`
2. Create project: `pr = Project("my_project")`
3. Create job: `job = pr.create_job(pr.job_type.Lammps, "job_name")`
4. Setup: `job.structure = ...; job.calc_minimize()`
5. Run: `job.run()`
6. Analyze: `job.output.energy_pot[-1]`

## Performance Characteristics
- Python overhead
- Efficient data I/O via HDF5
- Simplifies management of thousands of jobs

## Application Areas
- High-throughput thermodynamics (melting points, phase transitions)
- Defect energetics
- Machine learning potential training
- Interactive simulation tutorials

## Community and Support
- Developed by Max Planck Institut für Eisenforschung (MPIE)
- Active development
- Workshops and tutorials available

## Verification & Sources
**Primary sources**:
1. Homepage: https://pyiron.org/
2. GitHub: https://github.com/pyiron/pyiron
3. Publication: J. Janssen et al., Comp. Mater. Sci. 163, 24 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (MPIE)
- Applications: ICME, workflow IDE, Jupyter integration
