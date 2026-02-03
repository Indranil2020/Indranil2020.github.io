# GASP (Genetic Algorithm for Structure and Phase Prediction)

## Official Resources
- **Homepage**: http://gasp.mse.cornell.edu/
- **Source Repository**: https://github.com/henniggroup/GASP-python
- **Documentation**: https://github.com/henniggroup/GASP-python/blob/master/manual/manual.pdf
- **License**: GPL-3.0

## Overview
GASP (Genetic Algorithm for Structure and Phase Prediction) is a Python-based evolutionary algorithm package developed by the Hennig Group (Cornell/University of Florida). It is designed to predict stable crystal structures and phase diagrams by interfacing with ab initio (VASP) or classical (LAMMPS, GULP) energy calculators.

**Scientific domain**: Crystal structure prediction, phase diagram determination, evolutionary algorithms  
**Target user community**: Materials scientists, physicists, method developers

## Theoretical Methods
- **Genetic Algorithm (GA)**: Global optimization strategy.
- **Grand Canonical GA**: Allows variable composition searches for phase diagram prediction.
- **Evolutionary Operators**: Crossover, mutation, permutation, strain.
- **Energy Evaluation**: External interfaces to DFT (VASP) or classical potentials.

## Capabilities
- **Crystal Structure Prediction**: Finds low-energy structures for fixed compositions.
- **Phase Diagram Prediction**: Variable composition search to identify stable stoichiometries (convex hull).
- **Potentials Testing**: Can be used to fit or test empirical potentials against DFT data.
- **Interfacing**: Supports VASP, LAMMPS, GULP, and MOPAC.
- **Symmmetrization**: Can enforce or detect symmetry in generated structures.

## Inputs & Outputs
- **Input**:
  - `gasp_input.xml` or similar configuration file.
  - Calculator input files (e.g., `INCAR`, `POTCAR` for VASP).
- **Output**:
  - `run_data` directory containing structure files (POSCAR format).
  - `statistics` files tracking energy and evolution.
  - `best_structures` list.

## Interfaces & Ecosystem
- **VASP**: Primary first-principles engine.
- **LAMMPS/GULP**: For classical forcefield calculations.
- **Python**: Written in Python 2.7 (legacy) / Python 3 (modern branches).

## Verification & Sources
- **Confidence**: ✅ VERIFIED
- **Primary Source**: [GASP GitHub Repository](https://github.com/henniggroup/GASP-python)
- **Reference**: W. W. Tipton and R. G. Hennig, "A grand canonical genetic algorithm for the prediction of multi-component phase diagrams and testing of empirical potentials", *J. Phys.: Condens. Matter* 25, 495401 (2013).
