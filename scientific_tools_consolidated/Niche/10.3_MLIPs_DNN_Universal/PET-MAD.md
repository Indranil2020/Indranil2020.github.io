# PET-MAD

## Official Resources
- Source Repository: https://github.com/lab-cosmo/pet-mad
- License: Open source (MIT)

## Overview
**PET-MAD** is a lightweight universal interatomic potential trained on r²SCAN data. Covering 94 elements (102 in v1.5), it achieves competitive accuracy with only 2.8M parameters, making it efficient for production MD.

**Scientific domain**: Lightweight universal potential on r²SCAN data  
**Target user community**: Researchers needing efficient universal potential with meta-GGA accuracy

## Theoretical Methods
- Point Edge Transformer (PET) architecture
- MAD dataset (r²SCAN-based)
- Lightweight model (2.8M parameters)
- 94-102 elements coverage
- Metatrain fine-tuning

## Capabilities (CRITICAL)
- 94 elements (102 in v1.5)
- r²SCAN training data
- 2.8M parameters (lightweight)
- Metatrain fine-tuning
- ASE, i-PI, LAMMPS interfaces

**Sources**: GitHub repository, arXiv:2503.02089

## Key Strengths

### r²SCAN Training:
- Meta-GGA quality data
- Better than PBE-trained models
- Higher accuracy reference
- Chemical accuracy potential

### Lightweight:
- 2.8M parameters
- Efficient inference
- Fast MD
- Low memory

### Fine-Tuning:
- Metatrain framework
- Custom dataset adaptation
- Transfer learning
- Active learning

## Inputs & Outputs
- **Input formats**: Structures (ASE)
- **Output data types**: Energies, forces, stresses

## Interfaces & Ecosystem
- **ASE**: Calculator
- **i-PI**: MD interface
- **LAMMPS**: MD engine
- **Metatrain**: Fine-tuning

## Performance Characteristics
- **Speed**: Fast (lightweight)
- **Accuracy**: r²SCAN-level
- **System size**: 1-10000+ atoms

## Computational Cost
- **MD**: Very efficient
- **Fine-tuning**: Minutes to hours

## Limitations & Known Constraints
- **New project**: Still maturing
- **r²SCAN data**: Limited dataset size
- **Lightweight**: May sacrifice some accuracy

## Comparison with Other Codes
- **vs MACE-MP-0**: PET-MAD is r²SCAN, MACE is PBE
- **vs CHGNet**: PET-MAD is lighter, CHGNet has charge
- **Unique strength**: Lightweight universal potential (2.8M params) trained on r²SCAN meta-GGA data

## Application Areas

### Production MD:
- Efficient MD simulations
- Meta-GGA quality at PBE cost
- Structure relaxation
- Energy screening

### Fine-Tuning:
- Custom chemistry
- Active learning
- Multi-fidelity

## Best Practices
- Use for r²SCAN-quality MD
- Fine-tune with Metatrain
- Compare with PBE-trained models

## Community and Support
- Open source (MIT)
- Lab Cosmo maintained
- Metatensor ecosystem

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/lab-cosmo/pet-mad

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Lightweight universal potential (2.8M params) trained on r²SCAN meta-GGA data
