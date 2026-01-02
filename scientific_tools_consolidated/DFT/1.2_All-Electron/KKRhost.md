# KKRhost (Korringa-Kohn-Rostoker Host Code)

## Official Resources
- Homepage: Part of JuKKR suite - https://jukkr.fz-juelich.de/
- Documentation: https://jukkr.fz-juelich.de/
- Source Repository: https://iffgit.fz-juelich.de/kkr (Jülich GitLab)
- License: Academic/research (Forschungszentrum Jülich)

## Overview
KKRhost is the host system calculation module within the JuKKR suite, implementing the full-potential Korringa-Kohn-Rostoker (KKR) Green's function method for crystalline solids. Developed at Forschungszentrum Jülich, KKRhost calculates the electronic structure of perfect periodic crystals (the "host" system), providing the reference state for subsequent impurity calculations with KKRimp. It is part of the integrated JuKKR package for accurate all-electron DFT calculations.

**Scientific domain**: All-electron DFT, KKR Green's function method, crystalline solids  
**Target user community**: Materials scientists, condensed matter physicists, researchers studying periodic systems

## Theoretical Methods
- Korringa-Kohn-Rostoker (KKR) method
- Green's function formalism
- Multiple scattering theory
- All-electron approach
- Full-potential implementation
- Density Functional Theory
- LDA, GGA functionals
- Spin-polarized calculations
- Relativistic effects (scalar, full)

## Capabilities (CRITICAL)
**Category**: Academic/research code (Jülich)
- Ground-state electronic structure (perfect crystals)
- All-electron calculations
- Full-potential KKR
- Host system reference states
- Band structure calculations
- Density of states
- Magnetic properties
- Spin-polarized DFT
- Relativistic calculations
- Green's function methods
- Integration with KKRimp for impurities
- Part of JuKKR suite

**Sources**: JuKKR documentation and Jülich materials

## Key Strengths

### KKR Green's Function:
- Multiple scattering theory
- Efficient for complex geometries
- Natural for impurity problems
- Green's function advantages
- Spectroscopic applications

### All-Electron Accuracy:
- No pseudopotentials
- Full electronic structure
- Core state treatment
- Relativistic effects
- High accuracy

### Host-Impurity Framework:
- Perfect crystal host (KKRhost)
- Impurity calculations (KKRimp)
- Efficient defect studies
- Local perturbations
- Integrated workflow

### JuKKR Integration:
- Part of comprehensive suite
- Consistent methodology
- Shared infrastructure
- Multiple modules
- Unified framework

## Inputs & Outputs
- **Input formats**:
  - Crystal structure (conventional cells)
  - Atomic positions
  - KKR-specific parameters
  - Convergence settings
  
- **Output data types**:
  - Electronic structure
  - Green's functions
  - Band structure
  - Density of states
  - Magnetic moments
  - Host reference for KKRimp

## Interfaces & Ecosystem
- **JuKKR Suite**:
  - KKRhost (this code) - host systems
  - KKRimp - impurity calculations
  - KKRnano - large-scale KKR
  - Integrated workflow
  
- **Jülich Infrastructure**:
  - HPC support
  - JUDFT framework
  - AiiDA integration (potential)
  - European HPC access

## Workflow and Usage

### Host System Calculation:
```bash
# Calculate perfect crystal host
kkrhost < input.in > output.log

# Generate reference Green's functions
# Provide host state for impurity calculations
```

### Typical Workflow:
1. Define crystal structure (host)
2. Set KKR parameters
3. Run self-consistent calculation
4. Obtain host Green's functions
5. Use as reference for KKRimp defect studies

### Integration with KKRimp:
1. KKRhost: Calculate perfect crystal
2. Output: Host Green's functions
3. KKRimp: Add impurity/defect
4. Calculate local electronic structure changes

## Advanced Features

### Full-Potential KKR:
- Beyond muffin-tin
- Accurate charge density
- Complex materials
- High precision

### Relativistic Options:
- Scalar relativistic
- Fully relativistic
- Spin-orbit coupling
- Heavy elements

### Magnetic Systems:
- Spin-polarized DFT
- Magnetic moments
- Non-collinear magnetism
- Magnetic materials

### Green's Function:
- Energy-resolved
- Complex energy integration
- Spectral properties
- Efficient methods

## Performance Characteristics
- **Speed**: Efficient for KKR applications
- **Accuracy**: All-electron high accuracy
- **System size**: Moderate (periodic systems)
- **Purpose**: Reference host calculations
- **Typical**: Part of host-impurity workflow

## Computational Cost
- All-electron: More expensive than pseudopotential
- KKR method: Scales differently from plane-wave
- Efficient: For defect/impurity workflows
- HPC: Suitable for parallel computing

## Limitations & Known Constraints
- **Availability**: Academic access (Jülich)
- **Learning curve**: KKR methodology expertise
- **Documentation**: Research-level
- **Community**: Specialized (KKR users)
- **System types**: Primarily crystalline hosts
- **Distribution**: Limited compared to open-source codes

## Comparison with Other Codes
- **vs AkaiKKR**: JuKKR more modern, full-potential
- **vs SPR-KKR**: Similar methodology, different implementation
- **vs FLAPW codes**: KKR different formalism, advantages for defects
- **Unique strength**: Host-impurity framework, Green's function approach, JuKKR integration

## Application Areas

### Perfect Crystals:
- Host system calculations
- Reference electronic structure
- Bulk properties
- Periodic materials

### Defect Studies Preparation:
- Host state for impurities
- Reference for perturbations
- Defect formation energies
- Local electronic structure

### Magnetic Materials:
- Magnetic host systems
- Spin-polarized calculations
- Magnetic impurities (with KKRimp)

### Spectroscopy:
- Host spectral functions
- Reference states
- Green's function methods

## Best Practices

### Host Calculations:
- Converge host carefully
- High-quality Green's functions
- Proper energy resolution
- Validate with experiments

### Integration with KKRimp:
- Consistent parameters
- Same basis sets
- Proper energy grids
- Smooth workflow

### Jülich Resources:
- Access HPC facilities
- Consult documentation
- Contact support team
- Collaborate with developers

## Community and Support
- Forschungszentrum Jülich
- JuKKR user community
- Academic access
- HPC support (Jülich)
- Research collaborations
- European network

## Educational Resources
- JuKKR documentation
- Jülich tutorials
- KKR method literature
- Green's function theory
- Academic publications
- Training workshops (Jülich)

## Development
- Forschungszentrum Jülich
- JUDFT team
- Active development
- German research funding
- European collaborations
- HPC optimization

## Research Context
Part of comprehensive KKR framework at Jülich for accurate electronic structure calculations, particularly suited for impurity and defect studies where the host-impurity decomposition provides computational efficiency and physical insight.

## Verification & Sources
**Primary sources**:
1. JuKKR homepage: https://jukkr.fz-juelich.de/
2. Jülich GitLab: https://iffgit.fz-juelich.de/kkr
3. Forschungszentrum Jülich documentation
4. JUDFT team publications

**Secondary sources**:
1. KKR method literature
2. Green's function DFT papers
3. JuKKR publications
4. Jülich scientific reports

**Confidence**: VERIFIED - Academic code (Jülich)

**Verification status**: ✅ VERIFIED
- Institution: Forschungszentrum Jülich
- Access: Academic/research
- **Category**: Academic research code
- Status: Actively maintained
- Community: KKR/JuKKR users
- Specialized strength: KKR Green's function method for host systems, all-electron accuracy, integration with KKRimp for defects, part of comprehensive JuKKR suite, Jülich HPC infrastructure
