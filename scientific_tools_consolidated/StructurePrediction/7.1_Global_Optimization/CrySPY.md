# CrySPY

## Overview
CrySPY (pronounced "crispy") is a crystal structure prediction tool written in Python. It supports multiple search algorithms including random search (RS), Bayesian optimization (BO), Look Ahead based on Quadratic Approximation (LAQA), and evolutionary algorithms (EA).

## Theoretical Basis
- Random structure searching with symmetry constraints
- Bayesian optimization for efficient sampling
- LAQA for accelerated local optimization
- Evolutionary/genetic algorithms
- Interface with ML potentials for fast screening

## Key Capabilities
- Multiple CSP algorithms in one package
- Interface with VASP, QE, soiap, LAMMPS
- Support for ML potentials (M3GNet, MACE)
- Automatic structure generation with PyXtal
- Parallel execution support

**Sources**: CrySPY documentation, Sci. Technol. Adv. Mater. Methods 1, 87 (2021)

## Key Strengths

### Algorithm Variety:
- Random search, Bayesian optimization
- LAQA, evolutionary algorithms
- Flexible algorithm selection

### Interfaces:
- VASP, Quantum ESPRESSO
- LAMMPS, soiap
- ML potentials (M3GNet, MACE)

### Usability:
- Python-based, easy installation
- Good documentation
- Active development

## Inputs & Outputs
- **Input formats**: cryspy.in (configuration), initial structures (optional)
- **Output data types**: Optimized structures, energy rankings, logs

## Interfaces & Ecosystem
- **DFT codes**: VASP, Quantum ESPRESSO, OpenMX
- **ML potentials**: M3GNet, MACE, CHGNet
- **Structure generation**: PyXtal integration

## Workflow and Usage
1. Prepare cryspy.in configuration file
2. Set up calculator (VASP/QE/ML potential)
3. Run: `cryspy`
4. Monitor progress and collect results
5. Analyze lowest energy structures

## Performance Characteristics
- Efficient with ML potential pre-screening
- Scales with number of atoms and algorithm choice
- Bayesian optimization reduces required evaluations

## Computational Cost
- Depends on calculator (DFT vs ML)
- BO significantly reduces iterations
- Parallelizable across structures

## Best Practices
- Use ML potentials for initial screening
- Validate with DFT for final structures
- Choose algorithm based on system complexity
- Use symmetry constraints when appropriate

## Limitations & Known Constraints
- Performance depends on calculator choice
- Large systems require ML pre-screening
- Algorithm selection requires experience

## Application Areas
- Inorganic crystal structure prediction
- High-pressure phase discovery
- Materials screening with ML potentials
- Alloy structure prediction

## Comparison with Other Codes
- **vs USPEX**: CrySPY more algorithms, USPEX more mature
- **vs CALYPSO**: CrySPY Python-based, easier customization
- **vs XtalOpt**: CrySPY more ML integration
- **Unique strength**: Multiple algorithms, ML potential integration, Python flexibility

## Community and Support
- Open-source (MIT License)
- GitHub repository with active development
- Documentation and tutorials available
- Growing user community

## Verification & Sources
**Primary sources**:
1. Homepage: https://tomoki-yamashita.github.io/CrySPY_doc/
2. GitHub: https://github.com/Tomoki-YAMASHITA/CrySPY
3. Publication: Sci. Technol. Adv. Mater. Methods 1, 87 (2021)

**Secondary sources**:
1. CrySPY tutorials
2. PyXtal documentation
3. Published applications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub, MIT)
- Development: ACTIVE
- Applications: Crystal structure prediction, ML-accelerated CSP
