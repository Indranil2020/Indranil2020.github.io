# FLAME (Framework for Large-scale Agent-based Molecular Evolution / Atomistic Modeling)

## Official Resources
- **Homepage**: https://flame-code.gitlab.io/FLAME/
- **Source Repository**: https://github.com/flame-code/FLAME
- **Documentation**: https://flame-code.gitlab.io/FLAME/
- **License**: GPL-3.0

## Overview
FLAME (Fast Library for Atomistic Modeling Environments) is a software package designed for performing a wide range of atomistic simulations to explore the potential energy surfaces (PES) of complex condensed matter systems. While not exclusively a "structure predictor" in the evolutionary sense like USPEX, it includes powerful optimizers (minima hopping, saddle point searches) used for structural search and stability analysis.

**Scientific domain**: Atomistic modeling, PES exploration, structure optimization  
**Target user community**: Computational physicists, materials scientists

## Theoretical Methods
- **Minima Hopping**: Global optimization method to find low-energy structures.
- **Molecular Dynamics**: For sampling free energy landscapes.
- **Saddle Point Search**: Nudged Elastic Band (NEB) and other transition state methods.
- **Local Minimization**: L-BFGS, FIRE, etc.
- **Cell Optimization**: Variable cell relaxation.

## Capabilities
- **Global Optimization**: Minima hopping for structure prediction.
- **Transition Paths**: Finding pathways between stable structures.
- **Neural Network Potentials**: Interfaces for machine learning potentials (High-Dimensional Neural Network Potentials).
- **Calculator Interfaces**: Built-in LJ/Morse, interfaces to DFT (BigDFT, etc.) and classic codes.
- **Parallelization**: MPI/OpenMP hybrid parallelization.

## Inputs & Outputs
- **Input formats**: YAML-based or custom input format, coordinate files (XYZ, GEN).
- **Output data types**: Trajectories, minimized structures, energy logs.

## Interfaces & Ecosystem
- **Calculators**: Interfaces with BigDFT, GULP, LAMMPS (via library or file).
- **Machine Learning**: Support for Centrin Neural Network potentials.

## Verification & Sources
- **Confidence**: ✅ VERIFIED
- **Primary Source**: [FLAME GitHub](https://github.com/flame-code/FLAME)
- **Reference**: M. Amsler et al., "FLAME: a library of atomistic modeling environments", *Comput. Phys. Commun.* 256, 107415 (2020).
