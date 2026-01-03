# hiPhive

## Official Resources
- Homepage: https://hiphive.materialsmodeling.org/
- Documentation: https://hiphive.materialsmodeling.org/documentation.html
- Source Repository: https://gitlab.com/materials-modeling/hiphive
- License: MIT License

## Overview
hiPhive is a Python library for efficiently extracting high-order force constants from ab-initio molecular dynamics or systematic displacement calculations. Using advanced statistical methods including compressive sensing and sparse regression, hiPhive constructs force constant models for accurate prediction of phonon and thermal properties including anharmonic effects. The code integrates seamlessly with ASE and is designed for high-throughput phonon calculations.

**Scientific domain**: Anharmonic lattice dynamics, force constant extraction, thermal transport  
**Target user community**: Computational materials scientists, phonon researchers, thermal transport community

## Theoretical Methods
- Harmonic and anharmonic force constants (2nd, 3rd, 4th order)
- Compressive sensing algorithms
- Sparse regression methods (LASSO, elastic net)
- Least-squares fitting
- Model selection and validation
- Cross-validation techniques
- Symmetry constraints
- Cluster expansion for force constants
- Temperature-dependent effective potential (TDEP)

## Capabilities (CRITICAL)
- Extract force constants from ab-initio MD or displacement calculations
- 2nd, 3rd, and 4th order force constants
- Compressive sensing for minimal training data
- Automatic model selection and validation
- Phonon dispersion including anharmonicity
- Integration with phonopy for harmonic phonons
- Integration with phono3py for thermal conductivity
- Systematic displacement generation
- Force constant model validation
- Compatibility with multiple DFT codes via ASE
- Python API for custom workflows
- High-throughput capability

**Sources**: Official hiPhive documentation, Adv. Theory Simul. 2, 1800184 (2019)

## Key Strengths
- **Efficient extraction**: Compressive sensing reduces required calculations significantly
- **High-order IFCs**: Supports 2nd, 3rd, 4th order for strong anharmonicity
- **Python integration**: Native Python, integrates with ASE/phonopy/phono3py ecosystem
- **Validation tools**: Built-in cross-validation and model selection

## Inputs & Outputs
- **Input formats**:
  - ASE Atoms objects
  - Ab-initio MD trajectories (via ASE)
  - Systematic displacement structures
  - Forces from any DFT code (via ASE)
  
- **Output data types**:
  - Force constant matrices (2nd, 3rd, 4th order)
  - Model parameters and validation metrics
  - Phonopy-compatible force constants
  - phono3py-compatible 3rd order IFCs

## Interfaces & Ecosystem
- **ASE**: Native integration for structure handling
- **phonopy**: Export harmonic force constants
- **phono3py**: Export anharmonic force constants for thermal conductivity
- **DFT codes**: Any code supported by ASE (VASP, QE, GPAW, etc.)
- **Python**: Full Python API for scripting

## Workflow and Usage

### Typical Workflow:
```python
from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
from hiphive.fitting import Optimizer
from ase.io import read

# 1. Define cluster space (cutoffs, orders)
cs = ClusterSpace(prim_structure, cutoffs=[8.0, 5.0, 4.0])

# 2. Add training structures
sc = StructureContainer(cs)
for structure in training_structures:
    sc.add_structure(structure)

# 3. Train model with compressive sensing
opt = Optimizer(sc.get_fit_data(), train_size=0.8)
opt.train()

# 4. Extract force constants
fcp = ForceConstantPotential(cs, opt.parameters)
fcs = fcp.get_force_constants(supercell)

# 5. Export to phonopy/phono3py
fcp.write_to_phonopy(supercell, 'FORCE_CONSTANTS')
```

### MD-based Training:
```python
from ase.io import read

# Read AIMD trajectory
structures = read('aimd.traj', ':')

# Add to training set
for s in structures:
    sc.add_structure(s)
```

## Advanced Features
- **Compressive sensing**: Minimal training data via sparse optimization
- **Cross-validation**: Automatic model selection and overfitting prevention
- **High-order**: 3rd and 4th order anharmonic force constants
- **TDEP compatibility**: Temperature-dependent effective potential extraction
- **Model comparison**: Compare different cutoff/order combinations

## Performance Characteristics
- **Speed**: Python-based; moderate speed, focus on accuracy
- **Training data**: Compressive sensing reduces required data by ~50-70%
- **Scalability**: Handles large supercells efficiently
- **Memory**: Moderate; depends on cluster expansion size

## Computational Cost
- DFT calculations (MD or displacements) most expensive
- hiPhive fitting: Minutes to hours
- Model validation: Fast
- Overall: Efficient due to reduced training requirements

## Limitations & Known Constraints
- **Python overhead**: Slower than compiled codes for very large systems
- **Training data**: Still requires ab-initio calculations (reduced but necessary)
- **Cutoff selection**: Requires careful testing
- **Model validation**: Critical to avoid overfitting
- **Documentation**: Good but assumes familiarity with Python/ASE

## Comparison with Other Codes
- **vs ALAMODE**: hiPhive more Python-centric, ALAMODE more established
- **vs phono3py direct**: hiPhive uses compressive sensing, potentially fewer calculations
- **Unique strength**: Compressive sensing for efficient force constant extraction

## Application Areas
- **Anharmonic phonons**: Materials with strong anharmonicity
- **Thermal conductivity**: Input for phono3py calculations
- **High-throughput**: Automated force constant extraction
- **Strongly anharmonic systems**: Requires high-order force constants
- **Machine learning potentials**: Training data generation

## Best Practices
- Systematic convergence testing of cutoffs and orders
- Cross-validation to prevent overfitting
- Sufficient training data diversity
- Validate against test structures not in training
- Use symmetry constraints appropriately

## Community and Support
- Open-source (MIT license)
- GitLab repository
- Documentation website
- Active development
- Growing user community
- Integration with ASE ecosystem

## Educational Resources
- Comprehensive documentation
- Tutorial examples
- Publication with methodology
- Python API documentation
- Example notebooks

## Development
- Materials Modeling group, Chalmers University
- Active development
- Regular updates
- Community contributions welcome

## Research Impact
hiPhive enables efficient extraction of high-order force constants using compressive sensing, reducing computational cost for anharmonic phonon calculations and thermal transport studies.

## Verification & Sources
**Primary sources**:
1. Homepage: https://hiphive.materialsmodeling.org/
2. Documentation: https://hiphive.materialsmodeling.org/documentation.html
3. GitLab: https://gitlab.com/materials-modeling/hiphive
4. Publication: Adv. Theory Simul. 2, 1800184 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitLab, MIT)
- Development: ACTIVE (Chalmers)
- Applications: Anharmonic force constant extraction, compressive sensing, Python/ASE integration, phonopy/phono3py compatibility, high-throughput phonon calculations
