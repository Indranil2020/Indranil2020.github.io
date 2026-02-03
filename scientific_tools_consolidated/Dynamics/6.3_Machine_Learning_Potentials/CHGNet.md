# CHGNet

## Official Resources
- Homepage: https://chgnet.lbl.gov/
- Documentation: https://github.com/CederGroupHub/chgnet
- Source Repository: https://github.com/CederGroupHub/chgnet
- License: BSD-3-Clause

## Overview
CHGNet (Crystal Hamiltonian Graph neural Network) is a pretrained universal neural network potential for charge-informed atomistic modeling. It is trained on over 1.5 million structures from the Materials Project and explicitly accounts for magnetic moments and charge states.

**Scientific domain**: Universal ML potentials, materials science, charge-informed modeling  
**Target user community**: Materials scientists needing universal potentials with charge information

## Theoretical Methods
- Graph neural networks
- Charge-informed features
- Magnetic moment prediction
- Universal potential training
- Materials Project data

## Capabilities (CRITICAL)
- Pretrained universal potential
- Charge state prediction
- Magnetic moment prediction
- Oxidation state awareness
- ASE calculator
- Structure relaxation
- MD simulations

## Key Strengths

### Charge Information:
- Oxidation states
- Magnetic moments
- Charge transfer
- Redox reactions

### Universal Coverage:
- Most periodic table
- Materials Project training
- Ready to use

## Inputs & Outputs
- **Input formats**:
  - ASE Atoms
  - Pymatgen structures
  
- **Output data types**:
  - Energies
  - Forces
  - Stresses
  - Charges
  - Magnetic moments

## Interfaces & Ecosystem
- **ASE**: Calculator
- **Pymatgen**: Structure handling
- **Materials Project**: Training data

## Advanced Features
- **Charge prediction**: Oxidation states
- **Magnetism**: Magnetic moments
- **Universal**: Broad element coverage
- **Fine-tuning**: Transfer learning
- **Relaxation**: Structure optimization

## Performance Characteristics
- Fast inference
- GPU acceleration
- Good accuracy for materials
- Pretrained ready to use

## Computational Cost
- Pretrained: No training needed
- Inference: Fast
- Fine-tuning: Hours
- Overall: Very efficient

## Best Practices
- Use pretrained model first
- Validate for your system
- Fine-tune if needed
- Check charge predictions

## Limitations & Known Constraints
- Materials focus
- May need fine-tuning
- Charge accuracy varies
- Active development

## Application Areas
- Battery materials
- Catalysis
- Oxidation/reduction
- Magnetic materials
- Materials screening

## Comparison with Other Codes
- **vs M3GNet**: CHGNet charge-aware, M3GNet 3-body focus
- **vs MACE**: CHGNet charge/magnetism, MACE higher accuracy
- **vs DeepMD-kit**: CHGNet pretrained universal, DeepMD custom training
- **Unique strength**: Charge and magnetic moment prediction, Materials Project training

## Community and Support
- Active development (Berkeley/Ceder group)
- GitHub issues
- Good documentation
- Materials Project integration

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/CederGroupHub/chgnet
2. B. Deng et al., Nat. Mach. Intell. 5, 1031 (2023)

**Secondary sources**:
1. CHGNet tutorials
2. Materials Project documentation
3. Battery materials applications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, BSD-3)
- Published in Nature Machine Intelligence
- Academic citations: >200 (rapid growth)
- Active development: Ceder group
