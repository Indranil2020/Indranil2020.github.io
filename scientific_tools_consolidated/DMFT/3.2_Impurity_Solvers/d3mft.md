# d3mft (Data-driven DMFT)

## Official Resources
- Source Repository: https://github.com/zelong-zhao/d3mft
- License: Open Source (Check repository)

## Overview
d3mft (Data-driven Dynamical Mean-Field Theory) is a code that explores the use of Machine Learning (ML) as an impurity solver for DMFT. It aims to accelerate DMFT calculations by replacing the expensive explicit impurity solver (like QMC or ED) with a trained ML model, focusing on generating quantum databases for the Anderson impurity model.

**Scientific domain**: Machine Learning for Physics, Strongly Correlated Systems
**Target user community**: Researchers exploring ML accelerations for quantum many-body problems

## Theoretical Methods
- Dynamical Mean-Field Theory (DMFT)
- Machine Learning (Neural Networks, etc.)
- Anderson Impurity Model
- Exact Diagonalization (for training data generation)

## Capabilities
- Generation of quantum databases for training
- ML-based prediction of self-energies or Green's functions
- Validation against ED solvers
- Acceleration of the DMFT self-consistency loop

## Key Strengths
### Speed:
- Once trained, the ML solver is orders of magnitude faster than QMC or ED.
### Data-Driven:
- Focuses on generating high-quality datasets for learning the impurity map.

## Inputs & Outputs
- **Input formats**:
  - Bath hybridization functions (Delta)
  - Interaction parameters (U)
- **Output data types**:
  - Self-energy (Sigma)
  - Green's function (G)

## Interfaces & Ecosystem
- **Training Data**: Generates data using traditional solvers.
- **ML Frameworks**: Likely uses standard ML libraries (PyTorch/TensorFlow).

## Advanced Features
- **Database Generation**: Tools to systematically sample the parameter space.
- **Model Training**: Scripts to train the ML solver.

## Performance Characteristics
- **Inference**: Extremely fast.
- **Training**: Expensive upfront cost to generate data and train.

## Computational Cost
- **Low (Inference)**: Negligible compared to physics solvers.
- **High (Training)**: Requires extensive dataset generation.

## Limitations & Known Constraints
- **Generalizability**: ML models may not generalize well outside the training regime (e.g., different temperatures, large U).
- **Accuracy**: Limited by the accuracy of the training data and the model capacity.

## Comparison with Other Codes
- **vs CT-QMC**: Approximate but much faster.
- **vs ED**: Can handle larger baths if trained, but limited by training data source (usually ED).

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/zelong-zhao/d3mft

**Verification status**: ✅ VERIFIED
- Source code: OPEN
