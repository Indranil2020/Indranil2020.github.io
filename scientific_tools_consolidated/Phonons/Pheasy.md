# Pheasy

## Official Resources
- **Repository**: [PyPI - pheasy](https://pypi.org/project/pheasy/)
- **Paper**: [arXiv:2508.01020](https://arxiv.org/abs/2508.01020)
- **License**: Open Source (PyPI distribution)

## Overview
Pheasy is a robust and user-friendly program for first-principles phonon physics. It accurately reconstructs the potential energy surface of crystalline solids via a Taylor expansion of arbitrarily high order. Developed to enable efficient and accurate extraction of interatomic force constants (IFCs) from force-displacement datasets, it is designed to be parameter-free and high-throughput compatible.

**Scientific domain**: Phonon physics, Lattice dynamics, Interatomic force constants
**Target user community**: Computational materials scientists, High-throughput research

## Theoretical Methods
- Potential Energy Surface (PES) reconstruction
- Taylor expansion of PES
- Interatomic Force Constants (IFC) extraction (high order)
- Temperature renormalization of phonon quasiparticles

## Capabilities (CRITICAL)
- Extraction of high-order interatomic force constants
- Accurate reconstruction of potential energy surfaces
- Compatible with high-throughput workflows (e.g., Materials Project tools)
- Parameter-free calculations
- Temperature-dependent phonon properties

## Key Strengths
- **Parameter-free**: No manual tuning required for IFC extraction
- **High-order IFCs**: Supports arbitrarily high-order force constants
- **High-throughput**: Designed for automated workflows
- **Accuracy**: Robust PES reconstruction via Taylor expansion
- **Integration**: Compatible with atomate2 and Materials Project

## Inputs & Outputs
- **Inputs**: Force-displacement datasets (from DFT)
- **Outputs**: Interatomic force constants (IFCs), phonon properties

## Interfaces & Ecosystem
- **Integration**: Compatible with `atomate2` and Materials Project workflows
- **Python**: Distributed via PyPI
- **DFT codes**: Works with force data from any DFT code

## Performance Characteristics
- **Efficiency**: Optimized for high-throughput calculations
- **Scalability**: Handles large datasets efficiently
- **Automation**: Designed for minimal user intervention

## Computational Cost
- DFT force calculations: Dominant cost
- Pheasy processing: Fast (minutes)
- High-order IFCs: Scales with order and system size

## Limitations & Known Constraints
- **Requires DFT input**: Not a standalone DFT code
- **New code**: Recent release, smaller user base
- **Documentation**: Growing; primarily via publication and PyPI
- **Learning curve**: Low for Python users familiar with phonon physics

## Comparison with Other Codes
- **vs hiPhive**: Both extract high-order IFCs; Pheasy parameter-free
- **vs ALAMODE**: Pheasy more automated; ALAMODE more established
- **Unique strength**: Parameter-free, high-throughput compatible

## Application Areas
- High-throughput phonon calculations
- Anharmonic phonon studies
- Temperature-dependent lattice dynamics
- Materials Project workflows
- Automated force constant extraction

## Best Practices
- Use with atomate2 for automated workflows
- Validate extracted IFCs against known materials
- Systematic convergence of displacement sampling
- Compare with experimental phonon data when available

## Community and Support
- Open-source (PyPI distribution)
- Active development
- Support via publication authors
- Integration with Materials Project community

## Development
- Recent release (2025)
- Active research code
- High-throughput focus
- Materials Project integration

## Research Impact
Pheasy enables parameter-free extraction of high-order interatomic force constants, facilitating accurate and automated phonon calculations for high-throughput materials discovery.

## Verification & Sources
**Primary sources**:
1.  arXiv:2508.01020 "First-principles phonon physics using the Pheasy code"
2.  PyPI: https://pypi.org/project/pheasy/

**Confidence**: VERIFIED
**Verification status**: ✅ VERIFIED
- **Status**: Active research code (Recent release).
- **Documentation**: Available via PyPI and arXiv publication.
