# HTOCSP (High-Throughput Organic Crystal Structure Prediction)

## Official Resources
- **Homepage**: https://github.com/MaterSim/HTOCSP
- **Source Repository**: https://github.com/MaterSim/HTOCSP
- **Documentation**: https://github.com/MaterSim/HTOCSP/blob/main/README.md
- **License**: MIT License

## Overview
HTOCSP is a Python-based framework specifically designed for the automated, high-throughput prediction of organic crystal structures. It integrates various open-source tools to generate, optimize, and rank crystal structures of organic molecules, addressing the challenge of polymorphism in pharmaceutical and organic materials.

**Scientific domain**: Organic crystal structure prediction, high-throughput screening  
**Target user community**: Pharmaceutical researchers, organic chemists, computational materials scientists

## Theoretical Methods
- **Random Structure Generation**: Uses PyXtal for symmetry-compliant generation.
- **Force Field Optimization**: Uses CHARMM/Amber force fields via LAMMPS or GULP for initial screening.
- **DFT Optimization**: Refinement using DFT (e.g., VASP, Quantum ESPRESSO).
- **Clustering/Ranking**: Structure matching to identify unique polymorphs.

## Capabilities
- **Automated Workflow**: From SMILES/molecule to ranked crystal structures.
- **Symmetry Handling**: Generates structures in common organic space groups.
- **Force Field Assignment**: Automated parameterization for organic molecules.
- **Polymorph Screening**: Identifies low-energy packing arrangements.
- **Integration**: Wraps PyXtal, RDKit, and optimization engines.

## Inputs & Outputs
- **Input formats**: Molecular structure (SMILES, XYZ), search parameters.
- **Output data types**: Ranked crystal structures (CIF/POSCAR), energy tables.

## Interfaces & Ecosystem
- **PyXtal**: Core engine for structure generation.
- **RDKit**: Molecule handling and conformation generation.
- **LAMMPS/GULP**: Classical optimization.
- **VASP/QE**: Ab initio optimization.

## Verification & Sources
- **Confidence**: ✅ VERIFIED
- **Primary Source**: [HTOCSP GitHub](https://github.com/MaterSim/HTOCSP)
- **Reference**: Paper associated with MaterSim group (Recent 2024/2025 work).
