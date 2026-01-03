# MLIP (Machine Learning Interatomic Potentials)

## Official Resources
- Homepage: https://mlip.skoltech.ru/
- Documentation: https://gitlab.com/shapeev/mlip-2
- Source Repository: https://gitlab.com/shapeev/mlip-2
- License: MIT License

## Overview
MLIP is a software package for constructing moment tensor potentials (MTP). MTPs are a class of machine learning potentials that use polynomial invariants as descriptors. They are known for being computationally efficient (faster than neural networks) while maintaining accuracy comparable to GAP or NNPs. MLIP includes tools for active learning (generating training sets on the fly).

**Scientific domain**: Machine learning potentials, active learning  
**Target user community**: Metallurgists, MD users

## Capabilities (CRITICAL)
- **MTP**: Moment Tensor Potentials (polynomial basis).
- **Active Learning**: Algorithm (D-optimality) to select new configurations for training during MD.
- **LAMMPS**: Interface for MD.
- **Speed**: 10-100x faster than typical NNPs.

**Sources**: MLIP GitLab, Mach. Learn.: Sci. Technol. 1, 045022 (2020)

## Inputs & Outputs
- **Input formats**: CFG format (structures with forces/energies/stresses)
- **Output data types**: `.mtp` potential files

## Interfaces & Ecosystem
- **LAMMPS**: Pair style provided.
- **VASP/QE**: Interfaces for active learning loops.

## Workflow and Usage
1. Initial training set.
2. Train MTP: `mlp train init.mtp train.cfg > trained.mtp`
3. Run MD with active learning checks.
4. If extrapolation detected, run DFT on new structures and retrain.

## Performance Characteristics
- Highly efficient evaluation.
- Active learning minimizes the number of expensive DFT calculations needed.

## Application Areas
- Alloy phase diagrams
- Diffusion
- Crystal structure prediction

## Community and Support
- Developed by Shapeev Group (Skoltech)
- Active user base

## Verification & Sources
**Primary sources**:
1. GitLab: https://gitlab.com/shapeev/mlip-2
2. Publication: A. V. Shapeev, Multiscale Model. Simul. 14, 1153 (2016)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitLab)
- Development: ACTIVE
- Applications: MTP, active learning
