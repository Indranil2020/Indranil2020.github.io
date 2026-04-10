# DPNEGF

## Official Resources
- Source Repository: https://github.com/DeePTB-Lab/dpnegf
- Documentation: Included in repository
- License: Open source

## Overview
**DPNEGF** (DeePTB-NEGF) is a Python package that integrates the Deep Learning Tight-Binding (DeePTB) approach with the Non-Equilibrium Green's Function (NEGF) method, establishing an efficient quantum transport simulation framework with first-principles accuracy. It enables fast quantum transport calculations using ML-trained tight-binding Hamiltonians.

**Scientific domain**: ML-accelerated quantum transport, DeePTB-NEGF  
**Target user community**: Researchers needing fast quantum transport simulations with DFT accuracy using machine learning

## Theoretical Methods
- Non-equilibrium Green's function (NEGF)
- Deep Learning Tight-Binding (DeePTB)
- Slater-Koster parameterization
- LCAO Kohn-Sham Hamiltonian
- Machine learning Hamiltonian
- Quantum transport in open-boundary systems

## Capabilities (CRITICAL)
- Quantum transport with ML-trained Hamiltonians
- DeePTB-SK (Slater-Koster) transport
- DeePTB-E3 (LCAO) transport
- Open-boundary system simulation
- Transmission function calculation
- I-V characteristics
- First-principles accuracy at TB speed

**Sources**: GitHub repository, npj Comput. Mater.

## Key Strengths

### ML-Accelerated Transport:
- DFT accuracy at TB speed
- Orders of magnitude faster than DFT+NEGF
- Systematic improvement with training
- Generalizable to new structures

### DeePTB Integration:
- Well-established ML-TB framework
- Trained on DFT data
- Environment-corrected Hamiltonians
- Multiple basis options

### NEGF Framework:
- Open-boundary conditions
- Self-energy calculation
- Transmission and conductance
- Bias-dependent transport

## Inputs & Outputs
- **Input formats**:
  - DeePTB model files
  - Device structure
  - Transport configuration
  
- **Output data types**:
  - Transmission vs energy
  - I-V characteristics
  - Local density of states
  - Current density

## Interfaces & Ecosystem
- **DeePTB**: ML tight-binding framework
- **PyTorch**: ML backend
- **Python**: Scripting
- **NumPy**: Numerical computation

## Performance Characteristics
- **Speed**: Much faster than DFT+NEGF
- **Accuracy**: Near DFT quality
- **System size**: Thousands of atoms
- **Memory**: Moderate

## Computational Cost
- **ML prediction**: Fast (milliseconds per Hamiltonian)
- **NEGF calculation**: Minutes
- **Full I-V**: Hours
- **vs DFT+NEGF**: 10-100x speedup

## Limitations & Known Constraints
- **Training data dependent**: Quality limited by DeePTB training
- **DeePTB dependency**: Requires trained model
- **Newer code**: Less established than traditional NEGF
- **Documentation**: Limited

## Comparison with Other Codes
- **vs Transiesta**: DPNEGF is ML-fast, Transiesta is DFT-accurate but slow
- **vs Nanodcal**: DPNEGF is open source ML, Nanodcal is commercial LCAO
- **vs NanoNet**: DPNEGF uses ML Hamiltonians, NanoNet uses manual TB
- **Unique strength**: ML-accelerated quantum transport with DFT accuracy, DeePTB-NEGF integration

## Application Areas

### Nanoscale Electronics:
- FET device simulation
- Nanowire transport
- 2D material devices
- Contact resistance

### Break Junctions:
- Atomic-scale contacts
- Conductance quantization
- Structure-dependent transport
- Comparison with experiment

### High-Throughput Transport:
- Materials screening for devices
- Transport property databases
- Structure-property mapping
- Device optimization

## Best Practices

### DeePTB Training:
- Use diverse training structures
- Validate Hamiltonian accuracy
- Test extrapolation carefully
- Monitor energy convergence

### NEGF Calculation:
- Use sufficient energy grid
- Include adequate lead layers
- Check self-energy convergence
- Validate against DFT+NEGF

## Community and Support
- Open source on GitHub
- Developed by DeePTB Lab
- Published in npj Comput. Mater.
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/DeePTB-Lab/dpnegf
2. J. Zou et al., npj Comput. Mater. (2024)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Published methodology: npj Comput. Mater.
- Active development: Ongoing
- Specialized strength: ML-accelerated quantum transport with DFT accuracy, DeePTB-NEGF integration
