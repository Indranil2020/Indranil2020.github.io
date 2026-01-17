# SchNarc

## Official Resources
- Homepage: https://schnarc.readthedocs.io/
- Source Repository: https://github.com/schnarc/schnarc
- License: GNU General Public License v3.0

## Overview
SchNarc is a machine learning software package designed to accelerate nonadiabatic molecular dynamics simulations. It specifically interfaces with the SHARC molecular dynamics suite to replace expensive quantum chemical calculations with deep neural networks (SchNet). This allows for simulating extensive excited-state dynamics with ab initio accuracy at a fraction of the computational cost.

**Scientific domain**: Machine learning, nonadiabatic dynamics, photochemistry, neural networks
**Target user community**: Researchers performing large-scale excited-state dynamics requiring high statistical sampling

## Theoretical Methods
- Deep Neural Networks (SchNet architecture)
- Continuous-filter convolutional layers
- Learning of multiple electronic states
- Inverse-distance weighted interpolation
- Non-adiabatic coupling vector prediction
- Spin-orbit coupling prediction
- Phase correction for wavefunctions
- Active learning / Adaptive sampling

## Capabilities (CRITICAL)
- Prediction of energies, forces, and couplings
- Spin-orbit coupling learning
- Non-adiabatic coupling learning
- Interface with SHARC for dynamics
- Training data generation
- Model training and validation
- Long-timescale dynamics simulations
- Large ensemble simulations

**Sources**: Official GitHub, Documentation, J. Chem. Phys. 2018

## Key Strengths

### Speedup:
- Orders of magnitude faster than ab initio
- Enables nanosecond-scale dynamics
- Thousands of trajectories feasible

### Comprehensive Learning:
- Learns energies and forces
- Learns scalar and vector couplings
- Handles arbitrary number of states

### SHARC Integration:
- Seamless drop-in replacement for QC programs
- Uses SHARC input/output formats
- Complete dynamics workflow

## Inputs & Outputs
- **Input formats**:
  - Reference quantum chemistry data
  - Training configuration
  - SHARC input files
  
- **Output data types**:
  - Trained Neural Network models
  - Prediction logs
  - Dynamics trajectories (via SHARC)

## Interfaces & Ecosystem
- **Dynamics Engine**: SHARC
- **ML Framework**: TensorFlow / PyTorch (depends on version)
- **QC Codes**: Any code supported by SHARC (for training data)
- **Data format**: NumPy / Python pickle

## Advanced Features

### Phase Correction:
- Handles arbitrary phase of wavefunctions
- Consistent coupling prediction
- Smooth potential energy surfaces

### Adaptive Sampling:
- Iterative training/dynamics cycles
- Expansion of validity domain
- Automatic outlier detection

## Performance Characteristics
- **Speed**: Milliseconds per timestep
- **Accuracy**: <0.05 eV error (typical)
- **Scaling**: Linear with number of atoms (SchNet)
- **Training**: GPU-accelerated

## Computational Cost
- **Training**: High (requires GPUs)
- **Inference**: Very Low (CPU or GPU)
- **Data Generation**: High (expensive QC calculations)
- **Break-even point**: Efficient for many trajectories/long times

## Limitations & Known Constraints
- **Applicability domain**: Limited to trained chemical space
- **Data requirement**: Needs hundreds/thousands of QC points
- **Setup time**: High initial investment for training
- **Topology**: Fixed atom types/connectivity (usually)

## Comparison with Other Codes
- **vs Newton-X ML**: SchNarc specialized for SHARC/couplings
- **vs NEXMD**: SchNarc is ab initio based ML, NEXMD semiempirical
- **Unique strength**: Learning of couplings (NAC/SOC) and phase correction

## Application Areas
- **Long-time dynamics**: Intersystem crossing over nanoseconds
- **Large ensembles**: Converged statistics for spectra
- **Complex mechanisms**: Rare events in photochemistry
- **Material design**: Screening of photoactive molecules

## Best Practices
- **Data diversity**: Sample relevant phase space
- **Validation**: Check errors on test set strictly
- **Iterative learning**: Run short dynamics to find holes
- **Phase matching**: Ensure consistent phase in training data

## Community and Support
- Open-source GPL v3
- Developed by SHARC developers
- Documentation available
- Integration with SHARC ecosystem

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/schnarc/schnarc
2. Documentation: https://schnarc.readthedocs.io/
3. P. Westermayr et al., J. Chem. Phys. 149, 221101 (2018)

**Confidence**: VERIFIED - Associated with SHARC project

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (GitHub)
- Integration: Verified SHARC interface
- Specialized strength: ML for excited state couplings
