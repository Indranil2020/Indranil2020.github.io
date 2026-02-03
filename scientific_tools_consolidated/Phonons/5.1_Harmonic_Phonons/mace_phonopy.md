# mace_phonopy

## Official Resources
- Homepage: https://github.com/Mofahdi/mace_phonopy
- Source Repository: https://github.com/Mofahdi/mace_phonopy
- License: Open Source

## Overview
mace_phonopy is a code to generate second-order interatomic force constants (IFCs) from MACE machine learning potentials for use with Phonopy. It bridges MACE ML potentials with standard phonon calculation workflows.

**Scientific domain**: Machine learning potentials, phonon calculations, force constants  
**Target user community**: Researchers using MACE potentials for phonon calculations

## Theoretical Methods
- MACE machine learning potentials
- Finite displacement method
- Force constant extraction
- Phonopy integration
- Stability analysis

## Capabilities (CRITICAL)
- Generate FORCE_CONSTANTS from MACE
- Phonopy-compatible output
- Stability checking
- Phonon dispersion plotting
- band.conf generation
- Automated workflow

## Key Strengths

### MACE Integration:
- Direct MACE potential support
- Fast force evaluation
- Near-DFT accuracy
- GPU acceleration

### Phonopy Compatibility:
- Standard FORCE_CONSTANTS format
- band.conf generation
- Seamless integration
- Automated workflow

### Stability Analysis:
- Automatic stability check
- Imaginary frequency detection
- Stability criteria customization

## Inputs & Outputs
- **Input formats**:
  - Crystal structure (ASE-readable)
  - MACE potential file
  - Supercell parameters
  
- **Output data types**:
  - FORCE_CONSTANTS file
  - Stability status
  - band.conf file
  - Phonon dispersion plot

## Interfaces & Ecosystem
- **MACE**: ML potential
- **Phonopy**: Phonon calculations
- **ASE**: Structure handling
- **PyTorch**: ML backend

## Performance Characteristics
- **Speed**: Fast with GPU (MACE)
- **Accuracy**: Near-DFT quality
- **Scalability**: Large systems feasible

## Limitations & Known Constraints
- Requires trained MACE potential
- MACE-specific
- Potential quality dependent
- Requires PyTorch

## Comparison with Other Codes
- **vs Phonopy direct**: mace_phonopy uses ML forces; Phonopy uses DFT
- **vs autoplex**: Different ML potential (MACE vs MTP/GAP)
- **Unique strength**: Direct MACE-to-Phonopy bridge

## Best Practices

### Potential Quality:
- Use well-trained MACE potential
- Validate against DFT phonons
- Check force accuracy
- Test on known systems

### Calculations:
- Use appropriate supercell size
- Check convergence
- Validate stability results

## Application Areas
- High-throughput phonon screening
- Large system phonons
- ML potential validation
- Rapid phonon estimation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/Mofahdi/mace_phonopy

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Active development
