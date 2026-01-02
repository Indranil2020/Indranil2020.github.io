# ADF (Amsterdam Density Functional)

## Official Resources
- Homepage: https://www.scm.com/amsterdam-modeling-suite/adf/
- Documentation: https://www.scm.com/doc/
- Source Repository: Proprietary (commercial license)
- License: Commercial license (academic and commercial versions available)

## Overview
ADF (Amsterdam Density Functional) is a powerful DFT program particularly known for its Slater-type orbital (STO) basis sets, advanced relativistic methods, and spectroscopic property calculations. Developed by SCM (Software for Chemistry & Materials) in the Netherlands, ADF excels at molecular calculations, transition metal chemistry, spectroscopy, and accurate treatment of heavy elements. It is part of the Amsterdam Modeling Suite alongside BAND, DFTB, ReaxFF, and other modules.

**Scientific domain**: Molecular DFT, spectroscopy, relativistic quantum chemistry, transition metals  
**Target user community**: Chemists studying spectroscopy, catalysis, heavy elements, molecular properties

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Slater-type orbital (STO) basis sets
- Hybrid functionals (B3LYP, PBE0, etc.)
- Range-separated functionals
- Dispersion corrections (Grimme D3, D4)
- Time-Dependent DFT (TDDFT)
- Scalar relativistic methods (ZORA, X2C)
- Spin-orbit coupling
- Two-component and four-component relativistic
- Excited state gradients
- Conceptual DFT (Fukui functions, hardness)
- Fragment analysis and decomposition
- Solvation models (COSMO, SM12)

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and transition states
- Vibrational frequencies and IR/Raman spectra
- UV-Vis absorption and emission (TDDFT)
- Circular dichroism (CD, ECD, MCD)
- NMR chemical shifts and J-coupling
- EPR g-tensors and hyperfine coupling
- Mössbauer spectroscopy parameters
- X-ray absorption spectroscopy (XAS, XANES)
- Optical rotation and ORD
- Excited state dynamics
- Spin-spin coupling constants
- Electric field gradients (NQR)
- Polarizabilities and hyperpolarizabilities
- VCD (vibrational circular dichroism)
- Raman optical activity (ROA)
- Molecular orbitals and bonding analysis
- Energy decomposition analysis (EDA)
- Natural orbitals for chemical valence (NOCV)
- Atoms in molecules (AIM)
- Accurate heavy element chemistry
- Transition metal complexes
- Fragment-based calculations
- QM/MM methods

**Sources**: Official SCM documentation (https://www.scm.com/), confirmed in multiple source lists

## Key Strengths

### Slater-Type Orbitals:
- No basis set superposition error
- Accurate near nucleus
- Compact representation
- Natural for atoms
- Excellent for heavy elements

### Relativistic Methods:
- ZORA (zeroth-order regular approximation)
- Scalar and spin-orbit coupling
- X2C exact two-component
- Four-component Dirac
- Accurate for lanthanides/actinides

### Spectroscopy:
- Comprehensive spectroscopic properties
- NMR (chemical shifts, coupling)
- EPR (g-tensors, A-tensors)
- UV-Vis, CD, MCD
- XAS, Mössbauer
- VCD, ROA
- High accuracy

### Fragment Analysis:
- Energy decomposition analysis
- NOCV (deformation density)
- Charge transfer analysis
- Bonding understanding
- Conceptual DFT tools

### Transition Metals:
- Excellent for TM complexes
- Accurate d-orbital energies
- Spin states
- Ligand field effects
- Organometallics

## Inputs & Outputs
- **Input formats**:
  - Text-based input files
  - GUI (AMS-GUI) with visual builder
  - XYZ coordinates
  - PDB, MOL, SDF formats
  - Python scripting (PLAMS)
  
- **Output data types**:
  - Text output files
  - Binary data files
  - KF (Keyed File) format
  - Molecular orbitals
  - Densities and potentials
  - Spectra data
  - Formatted results

## Interfaces & Ecosystem
- **Amsterdam Modeling Suite**:
  - AMS-GUI (integrated interface)
  - BAND (periodic DFT)
  - DFTB (tight-binding)
  - ReaxFF (reactive force field)
  - Integrated workflow
  
- **Visualization**:
  - AMS-GUI built-in
  - ADFview
  - Export to standard formats
  - Molecular orbital visualization
  
- **Analysis**:
  - AMS analysis tools
  - Fragment analysis
  - Bonding analysis
  - Spectroscopy tools
  
- **Scripting**:
  - PLAMS (Python Library for Automating Molecular Simulation)
  - Python API
  - Workflow automation
  - High-throughput calculations

## Workflow and Usage

### GUI Workflow:
1. Build/import structure in AMS-GUI
2. Select ADF calculation type
3. Choose functional and basis set
4. Set calculation parameters
5. Run calculation
6. Visualize and analyze results

### Input File Example:

```
TITLE Water molecule

ATOMS
  O  0.0  0.0  0.0
  H  0.0  0.0  1.0
  H  0.0  1.0  0.0
END

BASIS
  Type TZP
END

XC
  GGA PBE
END

GEOMETRYOPTIMIZATION
END
```

### PLAMS Script:
```python
from scm.plams import *

init()
mol = Molecule('water.xyz')
sett = Settings()
sett.input.ams.Task = 'GeometryOptimization'
sett.input.adf.Basis.Type = 'TZP'
sett.input.adf.XC.GGA = 'PBE'
job = ADFJob(molecule=mol, settings=sett)
result = job.run()
finish()
```

## Advanced Features

### Energy Decomposition Analysis:
- Pauli repulsion
- Electrostatic interaction
- Orbital interactions
- Bonding understanding
- Fragment-based interpretation

### NOCV Analysis:
- Natural orbitals for chemical valence
- Deformation density
- Charge transfer channels
- σ/π bonding decomposition
- Visual interpretation

### Excited States:
- TDDFT for absorption/emission
- Excited state geometry optimization
- Radiative and non-radiative decay
- Phosphorescence
- Spin-orbit coupling effects

### Relativistic Calculations:
- ZORA (efficient, accurate)
- Spin-orbit coupling included
- Heavy element spectra
- Lanthanide/actinide chemistry
- Accurate bond energies

### Conceptual DFT:
- Fukui functions
- Chemical hardness/softness
- Electrophilicity
- Dual descriptors
- Reactivity indices

## Performance Characteristics
- **Speed**: Competitive for molecular systems
- **Accuracy**: Excellent for spectroscopy
- **System size**: Up to ~500 atoms practical
- **Memory**: Moderate requirements
- **Parallelization**: Good multi-core scaling

## Computational Cost
- **DFT**: Standard scaling
- **Hybrids**: More expensive
- **TDDFT**: Moderate cost
- **Relativistic**: Manageable overhead
- **Large molecules**: Feasible with modern hardware

## Limitations & Known Constraints
- **Commercial**: License required
- **Periodic systems**: Use BAND module instead
- **Very large systems**: Limited vs plane-wave codes
- **Cost**: Commercial licensing
- **STO basis**: Limited availability compared to Gaussians
- **Learning curve**: Moderate
- **Platform**: Windows, Linux, macOS

## Comparison with Other Codes
- **vs Gaussian**: ADF better for spectroscopy, heavy elements, STO basis
- **vs ORCA**: Both strong for spectroscopy, different approaches
- **vs NWChem**: ADF more user-friendly, better GUI
- **vs Turbomole**: Similar capabilities, different basis sets
- **Unique strength**: STO basis, comprehensive spectroscopy, fragment analysis, relativistic methods

## Application Areas

### Spectroscopy:
- NMR predictions
- EPR simulations
- UV-Vis spectra
- CD/MCD calculations
- XAS/XANES
- Vibrational spectroscopy

### Catalysis:
- Transition metal catalysts
- Reaction mechanisms
- Ligand effects
- Activation barriers
- Organometallic chemistry

### Heavy Elements:
- Lanthanides/actinides
- Relativistic effects
- Bonding in f-element compounds
- Nuclear properties

### Materials Chemistry:
- Molecular materials
- Optical properties
- Electronic structure
- Excited states

## Best Practices

### Basis Set Selection:
- DZ for quick tests
- DZP for standard
- TZP for publication quality
- TZ2P for high accuracy
- QZ4P for benchmark

### Functional Choice:
- LDA for quick tests
- GGA (PBE, BP86) for general use
- Hybrids (B3LYP, PBE0) for accuracy
- Range-separated for charge transfer
- Include dispersion for weak interactions

### Relativistic:
- Scalar ZORA for elements Z>36
- Spin-orbit for heavy elements
- X2C for highest accuracy
- Check frozen core approximation

### Convergence:
- Use good initial geometry
- Appropriate integration accuracy
- SCF convergence criteria
- Symmetry when applicable

### Spectroscopy:
- Include solvent effects
- Use appropriate functional
- Sufficient basis set
- Relativistic for heavy elements

## Community and Support
- Commercial support from SCM
- Comprehensive documentation
- Regular updates and new features
- Training courses available
- User community
- Email support

## Educational Resources
- Extensive documentation
- Tutorial examples
- Video tutorials
- Workshops and courses
- Application notes
- Published papers

## Development
- Active development by SCM
- Regular releases (annual)
- New features added
- Bug fixes and improvements
- User feedback incorporated
- Modern software architecture

## Amsterdam Modeling Suite Integration
- Seamless workflow between modules
- ADF for molecules
- BAND for periodic
- DFTB for speed
- ReaxFF for dynamics
- Unified interface

## Verification & Sources
**Primary sources**:
1. Official website: https://www.scm.com/amsterdam-modeling-suite/adf/
2. Documentation: https://www.scm.com/doc/
3. G. te Velde et al., J. Comput. Chem. 22, 931 (2001) - ADF methodology
4. E. van Lenthe et al., J. Chem. Phys. 99, 4597 (1993) - ZORA relativistic method

**Secondary sources**:
1. SCM documentation and tutorials
2. Published studies using ADF (>30,000 citations)
3. Spectroscopy validation papers
4. Confirmed in multiple source lists

**Confidence**: VERIFIED - Appears in multiple independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Software: Commercial (widely available)
- Community support: Excellent (SCM support, documentation)
- Academic citations: >35,000
- Active development: Regular annual releases
- Specialized strength: Slater-type orbitals, comprehensive spectroscopy, relativistic methods, fragment analysis, transition metal chemistry
