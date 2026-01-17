# Psi4NumPy

## Official Resources
- Homepage: https://github.com/psi4/psi4numpy
- Documentation: https://github.com/psi4/psi4numpy/tree/master/Tutorials
- Source Repository: https://github.com/psi4/psi4numpy
- License: BSD 3-Clause "New" or "Revised" License

## Overview
Psi4NumPy is an educational and development framework that bridges the Psi4 quantum chemistry package with the NumPy Python library. It provides interactive tutorials and reference implementations of modern quantum chemical methods. By exposing the core C++ modules of Psi4 to Python, it allows users to write clear, readable, and efficient implementations of complex methods like HF, MP2, CC, and CI directly in Python.

**Scientific domain**: Quantum Chemistry Education, Method Development  
**Target user community**: Students, Educators, Method Developers

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Møller-Plesset Perturbation Theory (MP2)
- Coupled Cluster (CCSD, CCSD(T))
- Configuration Interaction (CI)
- Symmetry-Adapted Perturbation Theory (SAPT)
- Density Functional Theory (DFT grids)
- Response Theory
- Geometry Optimization

## Capabilities (CRITICAL)
- Interactive Jupyter notebook tutorials
- Reference implementations of QC methods
- Access to Psi4 internals (Integrals, Wavefunctions)
- Matrix manipulation via NumPy
- Tensor contraction via opt_einsum
- Custom SCF solvers
- Property calculations
- educational visualization of algorithms

## Key Strengths

### Education:
- "Executable papers"
- Step-by-step algorithms
- Readable code vs optimized black-box
- Visualizing convergence and matrices
- Low barrier to entry

### Prototyping:
- Rapid method development
- Python flexibility with C++ backend
- Easy debugging
- Validation against production code

### Interoperability:
- NumPy ecosystem
- SciPy tools
- Matplotlib visualization
- Tensor libraries

## Inputs & Outputs
- **Input**: Jupyter Notebooks, Python scripts
- **Output**: Python objects, Arrays, Energies
- **Data**: Access to raw integrals and tensors

## Interfaces & Ecosystem
- **Psi4**: The core engine
- **NumPy**: The math engine
- **Jupyter**: The interface
- **Binder**: Cloud execution

## Advanced Features

### Tutorials:
- HF/DFT basics
- Response theory
- CEPA/CC theory
- Intermolecular forces (SAPT)
- Orbital rotations

### Developer Tools:
- Access to JK build objects
- DIIS implementations
- Orbital spaces
- Integral transformations

## Performance Characteristics
- **Speed**: Python overhead (slower than core Psi4)
- **Purpose**: Clarity over raw speed
- **Operations**: Heavy lifting by Psi4 core/BLAS
- **Scalability**: For small systems/learning

## Computational Cost
- **Learning**: Free (Open Source)
- **Compute**: Low (Small molecules)
- **Development**: Fast prototyping

## Limitations & Known Constraints
- **Performance**: Not for production runs on large systems
- **Scope**: Focuses on single-node educational/dev tasks
- **Dependency**: Requires Psi4

## Comparison with Other Codes
- **vs PySCF**: Similar Python focus, Psi4NumPy more tutorial-oriented
- **vs Standard Psi4**: Psi4NumPy is the *interface/tutorial* layer
- **vs Q-Chem**: Open and scriptable vs compiled commercial
- **Unique strength**: Best-in-class educational tutorials for QC coding

## Application Areas
### Academic Instruction:
- **Undergraduate QC**: Interactive demonstrations of orbitals and bonding
- **Graduate Courses**: Programming assignments for HF, MP2, and CI
- **Summer Schools**: Hands-on workshops for theory development
- **Self-Paced Learning**: Comprehensive set of progressively difficult tutorials

### Method Development:
- **Algorithm Prototyping**: Testing new density functionals or correlation methods
- **Tensor Operations**: Developing tensor contraction logic before C++ implementation
- **Validation**: Generating reference values for debugging compiled codes
- **Visualization**: Plotting convergence metrics and wavefunction properties interactively

## Best Practices
### Learning Path:
- **Beginner**: Start with `1_Psi4NumPy-Basics` to understand the data structures
- **Intermediate**: Proceed to `3_Hartree-Fock` to implement your own SCF code
- **Advanced**: Tackle `5_CEPA-and-CC` for correlated method internals

### Performance Tuning:
- **Vectorization**: Use NumPy broadcasting instead of Python loops where possible
- **Tensor Contraction**: Utilize `opt_einsum` for efficient tensor operations
- **Memory**: Be mindful of storing large rank-4 tensors (ERIs) in memory
- **Validation**: Check your implementation against Psi4's built-in energy values

## Community and Support
- Large GitHub community
- Psi4 developers
- Crawford group (Virginia Tech)
- Sherrill group (Georgia Tech)

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/psi4/psi4numpy
2. Smith et al., J. Chem. Theory Comput. 14, 3504 (2018)
3. Psi4Education methodology

**Confidence**: VERIFIED
- Status: Active open project
- Impact: widely used in education
