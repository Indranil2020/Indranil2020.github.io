# SchNetPack

## Official Resources
- Homepage: https://schnetpack.readthedocs.io/
- Documentation: https://schnetpack.readthedocs.io/
- Source Repository: https://github.com/atomistic-machine-learning/schnetpack
- License: MIT License

## Overview
SchNetPack is a deep learning toolbox for atomistic systems. It implements SchNet, a continuous-filter convolutional neural network, and other models (PaiNN) for modeling quantum chemical properties. It is designed to train potentials on large datasets of molecules and materials to predict energies, forces, and dipole moments with DFT accuracy but at a fraction of the cost.

**Scientific domain**: Machine learning potentials, deep learning, quantum chemistry  
**Target user community**: Computational chemists, ML researchers

## Capabilities (CRITICAL)
- **Models**: SchNet, PaiNN (equivariant networks).
- **Training**: Flexible training loops using PyTorch.
- **ASE Interface**: Can be used as an ASE calculator for MD.
- **Datasets**: Built-in access to MD17, QM9, Materials Project datasets.
- **Properties**: Energy, forces, dipole moment, polarizability, etc.

**Sources**: SchNetPack documentation, J. Chem. Theory Comput. 15, 448 (2019)

## Inputs & Outputs
- **Input formats**: ASE Atoms, XYZ, database files (ASE db)
- **Output data types**: Trained models (.pth), predicted properties

## Interfaces & Ecosystem
- **PyTorch**: Built on PyTorch.
- **ASE**: Integration for simulation.
- **PyTorch Lightning**: Used for training infrastructure (in v2.0+).

## Workflow and Usage
1. Load dataset (e.g., QM9).
2. Configure model (SchNet) and task (Energy/Forces).
3. Train model: `trainer.fit(model, datamodule)`.
4. Use for inference or MD.

## Performance Characteristics
- GPU-accelerated training and inference.
- Scales to large datasets (millions of conformations).
- State-of-the-art accuracy for organic molecules and bulk materials.

## Application Areas
- Molecular dynamics of large systems.
- Drug discovery (property prediction).
- Exploration of chemical space.

## Community and Support
- Developed by TU Berlin (Machine Learning Group)
- Active GitHub repository

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/atomistic-machine-learning/schnetpack
2. Publication: K. T. Schütt et al., J. Chem. Theory Comput. 15, 448 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Deep learning potentials, SchNet
