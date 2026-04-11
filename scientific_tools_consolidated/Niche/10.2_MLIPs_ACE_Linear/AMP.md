# AMP (Atomistic Machine-learning Package)

## Official Resources
- Homepage: https://amp.readthedocs.io/
- Documentation: https://amp.readthedocs.io/
- Source Repository: https://bitbucket.org/andrewpeterson/amp (migrated/archived)
- License: GPL v3

## Overview
AMP (Atomistic Machine-learning Package) is an open-source Python package for generating atomistic machine learning potentials (interatomic potentials). It interfaces with ASE to train potentials on DFT data (energy and forces) using descriptors like Gaussian symmetry functions and neural networks.

**Scientific domain**: Machine learning potentials, molecular dynamics  
**Target user community**: Catalysis researchers, MD users

## Capabilities (CRITICAL)
- **Training**: Train neural network potentials on DFT trajectories.
- **Descriptors**: Fingerprint schemes (Gaussian, Zernike).
- **ASE Interface**: Functions as an ASE calculator for MD simulations.
- **Regression**: Neural networks, Gaussian processes.

**Sources**: AMP documentation, Comp. Phys. Comm. 203, 154 (2016)

## Inputs & Outputs
- **Input formats**: ASE Trajectory files (training data)
- **Output data types**: Potential file (.amp)

## Interfaces & Ecosystem
- **ASE**: Deeply integrated.
- **TensorFlow/Scikit-learn**: Used for backend (legacy versions).

## Workflow and Usage
1. Generate training data (DFT MD).
2. `calc = Amp(descriptor=..., model=...)`
3. `calc.train(images=training_data)`
4. `atoms.calc = calc; atoms.get_forces()`

## Performance Characteristics
- Training is computationally intensive.
- Prediction is orders of magnitude faster than DFT.
- Speed comparable to other NN potentials.

## Application Areas
- Catalysis (surface reactions)
- Molecular dynamics of nanoparticles
- Global optimization

## Community and Support
- Developed by Peterson Group (Brown University)
- **Status**: Mostly legacy/stable; newer tools (SchNetPack, NequIP) are more active.

## Verification & Sources
**Primary sources**:
1. Documentation: https://amp.readthedocs.io/
2. Publication: A. A. Khorshidi and A. A. Peterson, Comp. Phys. Comm. 203, 154 (2016)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN
- Development: MAINTENANCE
- Applications: ML potentials, neural networks
