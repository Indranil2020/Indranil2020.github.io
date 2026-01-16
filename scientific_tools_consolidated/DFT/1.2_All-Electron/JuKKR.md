# JuKKR (Jülich KKR)

## Official Resources
- Homepage: https://jukkr.fz-juelich.de/
- Documentation: https://jukkr.readthedocs.io/
- Source Repository: https://github.com/JuDFTteam/JuKKR
- License: MIT License

## Overview
JuKKR is a modern KKR (Korringa-Kohn-Rostoker) Green's function DFT code developed at Forschungszentrum Jülich, Germany. Part of the JuDFT family of codes alongside Fleur, JuKKR provides comprehensive KKR method capabilities with modern Python interfaces (masci-tools, aiida-kkr) for high-throughput workflows and materials informatics. It implements the full-potential KKR method using multiple scattering theory for electronic structure, spectroscopy, and disordered alloys.

**Scientific domain**: All-electron KKR Green's function DFT, alloys, magnetism, spectroscopy  
**Target user community**: Materials scientists, alloy researchers, workflow developers, HPC users

## Theoretical Methods
- Korringa-Kohn-Rostoker (KKR) Green's function method
- Full-potential implementation
- Multiple scattering theory
- All-electron approach (no pseudopotentials)
- Coherent Potential Approximation (CPA) for disorder
- Density Functional Theory (LDA, GGA, meta-GGA)
- Scalar and fully relativistic treatments
- Spin-polarized and non-collinear magnetism
- Spin-orbit coupling

## Capabilities (CRITICAL)
- Ground-state electronic structure (all-electron)
- Full-potential KKR calculations
- Perfect crystalline solids (KKRhost module)
- Impurity and defect calculations (KKRimp module)
- Disordered alloys via Coherent Potential Approximation (CPA)
- Random substitutional alloys
- High-entropy alloys
- Band structure and density of states
- Magnetic properties (moments, anisotropies, exchange interactions)
- Spectroscopy calculations (XAS, XMCD, photoelectron spectroscopy)
- Transport properties
- Spin dynamics
- Green's function formalism for efficient calculations
- High-throughput screening (via AiiDA-KKR)
- Automated workflows and materials informatics
- Massively parallel calculations (KKRnano for extreme scale)

**Sources**: Official JuKKR website (https://jukkr.fz-juelich.de/), GitHub, Jülich documentation

## JuKKR Suite Components

### KKRhost
- Host system calculations for perfect periodic crystals
- Reference states for impurity calculations
- Full-potential all-electron DFT
- Magnetic and non-magnetic systems

### KKRimp
- Impurity and defect calculations
- Local perturbations in host materials
- Efficient Green's function embedding
- Defect formation energies

### KKRnano
- Massively parallel variant for extreme-scale HPC
- Hundreds of thousands of cores
- Large-scale materials simulations
- Leadership computing applications

## Key Strengths

### Modern Infrastructure
- Python-based workflows and interfaces
- AiiDA integration for reproducibility
- Automated high-throughput calculations
- Materials informatics ready
- Open-source (MIT License)

### KKR Green's Function Method
- Multiple scattering theory
- Natural for complex geometries
- Efficient for impurities and defects
- Direct access to spectroscopic properties
- No supercells needed for disorder (CPA)

### CPA for Alloys
- Random substitutional alloys
- Exact statistical treatment
- Concentration-dependent properties
- High-entropy alloys
- No supercell approximation

### JuDFT Ecosystem
- Part of comprehensive Jülich suite
- Integration with Fleur (FLAPW code)
- Shared tools and infrastructure
- Consistent methodology
- Strong HPC support

## Inputs & Outputs
- **Input formats**:
  - Python-based input generation (masci-tools)
  - Text-based input files
  - Crystal structure definitions
  - AiiDA workflow specifications
  
- **Output data types**:
  - Electronic structure data
  - Green's functions
  - Band structure and DOS
  - Magnetic properties
  - Spectroscopic data
  - AiiDA data provenance

## Interfaces & Ecosystem
- **Workflow Automation**:
  - AiiDA-KKR - High-throughput workflows
  - masci-tools - Python utilities for JuKKR
  - Automated convergence testing
  - Materials database integration
  
- **JuDFT Family Integration**:
  - Fleur (FLAPW) - complementary all-electron code
  - Shared Jülich infrastructure
  - Combined workflows
  
- **HPC Infrastructure**:
  - Forschungszentrum Jülich supercomputers
  - European HPC centers
  - Optimized for parallel computing

## Workflow and Usage

### Basic KKRhost Calculation
```bash
# Calculate perfect crystal host
kkrhost < input.in > output.log
```

### AiiDA-KKR Workflow (Python)
```python
from aiida import load_profile
from aiida_kkr.workflows import kkr_scf_wc

# Load AiiDA profile
load_profile()

# Set up structure and parameters
structure = ...  # AiiDA StructureData
parameters = {...}  # Calculation parameters

# Run self-consistent calculation
result = submit(kkr_scf_wc,
                structure=structure,
                parameters=parameters)
```

### CPA for Disordered Alloys
```bash
# Random alloy A_{x}B_{1-x}
# CPA treats disorder exactly without supercells
kkr_cpa < alloy_input.in > cpa_output.log
```

## Advanced Features

### Coherent Potential Approximation
- Substitutional disorder treatment
- Multiple components per site
- Concentration-dependent properties
- No supercell needed
- Exact statistical averaging

### Spectroscopy Capabilities
- X-ray absorption spectroscopy (XAS)
- X-ray magnetic circular dichroism (XMCD)
- Photoelectron spectroscopy
- Core-level excitations
- Element-specific magnetism

### Magnetic Properties
- Collinear and non-collinear magnetism
- Magnetic anisotropy calculations
- Exchange interactions (J_ij)
- Spin dynamics parameters
- Complex magnetic structures

## Performance Characteristics
- **Parallelization**: MPI and OpenMP support
- **Scalability**: Good for KKR applications, excellent with KKRnano
- **Typical systems**: Unit cells to moderate supercells
- **CPA efficiency**: No supercell overhead for disorder
- **HPC ready**: Optimized for Jülich supercomputers

## Limitations & Known Constraints
- **Learning curve**: KKR methodology requires understanding
- **All-electron cost**: More expensive than pseudopotential methods
- **Green's function complexity**: Different from plane-wave approaches
- **Documentation**: Comprehensive but technical
- **Best for**: Periodic systems, alloys, defects, spectroscopy

## Comparison with Other KKR Codes
- **vs AkaiKKR**: JuKKR more modern, full-potential, Python workflows
- **vs SPR-KKR**: Similar capabilities, JuKKR has AI iDA integration
- **vs Plane-wave DFT**: KKR natural for disorder/defects, different formalism
- **Unique strength**: Modern Python workflows, AiiDA integration, full JuDFT ecosystem, HPC-optimized

## Application Areas
- Disordered alloys and high-entropy alloys
- Magnetic materials and spintronics
- Defect and impurity physics
- Spectroscopy interpretation
- High-throughput materials screening
- Materials informatics and databases

## Community and Support
- Open-source (MIT License)
- Forschungszentrum Jülich support
- Active GitHub development
- AiiDA community integration
- European HPC network
- Workshops and tutorials

## Educational Resources
- Official documentation: https://jukkr.readthedocs.io/
- AiiDA-KKR tutorials
- KKR method literature
- Jülich training workshops
- GitHub examples and issues

## Development
- Forschungszentrum Jülich (IAS-1, PGI-1)
- JuDFT team
- Active open-source development
- Regular releases and updates
- Community contributions welcome

## Verification & Sources
**Primary sources**:
1. Official website: https://jukkr.fz-juelich.de/
2. GitHub: https://github.com/JuDFTteam/JuKKR
3. Documentation: https://jukkr.readthedocs.io/
4. AiiDA-KKR: https://github.com/JuDFTteam/aiida-kkr

**Secondary sources**:
1. Forschungszentrum Jülich publications
2. KKR method literature
3. JuDFT family documentation
4. AiiDA plugin registry

**Confidence**: VERIFIED - Official Jülich code

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- GitHub repository: OPEN (MIT License)
- Documentation: COMPREHENSIVE
- Active development: Regular updates
- Community support: GitHub, AiiDA forum, Jülich
- Specialized strength: Modern Python-based KKR with AiiDA workflows, full-potential Green's function method, CPA for alloys, HPC-optimized, comprehensive JuDFT ecosystem integration
