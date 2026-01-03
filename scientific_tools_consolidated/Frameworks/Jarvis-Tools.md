# Jarvis-Tools

## Official Resources
- Homepage: https://jarvis-tools.readthedocs.io/
- Documentation: https://jarvis-tools.readthedocs.io/
- Source Repository: https://github.com/usnistgov/jarvis
- License: NIST License

## Overview
jarvis-tools is the open-source Python software package that powers the JARVIS database. It provides a suite of tools for designing, managing, and analyzing atomistic simulations (DFT, MD) and applying machine learning to materials data. It supports VASP, Quantum ESPRESSO, LAMMPS, and other codes.

**Scientific domain**: Materials informatics, automation, machine learning  
**Target user community**: Users of JARVIS database, high-throughput researchers

## Capabilities (CRITICAL)
- **Automation**: Workflows for VASP (TB-mBJ, SOC), QE, and LAMMPS.
- **Analysis**: Band structure, DOS, topological invariants, elastic tensors, STM images.
- **Machine Learning**: Graph Convolutional Networks (GCN), descriptors, and pre-trained models.
- **Database**: Tools to interact with and download JARVIS datasets.
- **Wannier**: Tight-binding hamiltonian generation via Wannier90.

**Sources**: jarvis-tools documentation

## Inputs & Outputs
- **Input formats**: Atomic structures (POSCAR, CIF), calculation inputs
- **Output data types**: JSON, XML, Plots

## Interfaces & Ecosystem
- **Codes**: VASP, QE, LAMMPS, Wannier90, Boltztrap.
- **ML**: DGL, PyTorch, Scikit-learn.
- **JARVIS**: Official API.

## Workflow and Usage
1. Install: `pip install jarvis-tools`
2. Download data:
   ```python
   from jarvis.db.figshare import data
   d = data('dft_3d')
   ```
3. Run ML model: `AL = GraphConv(config).train(...)`

## Performance Characteristics
- Comprehensive suite for 2D and 3D materials
- Strong focus on post-processing and ML

## Application Areas
- High-throughput screening
- Machine learning model development
- Electronic structure analysis

## Community and Support
- Developed by NIST
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/usnistgov/jarvis

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: JARVIS ecosystem, ML, DFT automation
