# QuantumATK

## Official Resources
- Homepage: https://www.synopsys.com/silicon/quantumatk.html
- Documentation: https://docs.quantumatk.com/
- Classification: Commercial Software (Synopsys)

## Overview
QuantumATK (formerly Atomistix ToolKit) is a comprehensive software platform for atomic-scale modeling of materials, nanostructures, and devices. While it is famous for its NEGF transport capabilities using LCAO, it includes a robust **PlaneWave** calculator (PW) that allows for systematically improvable basis set calculations, making it a "complete" electronic structure suite.

**Scientific domain**: Materials science, semiconductor devices, quantum transport, electronic structure
**Target user community**: Semiconductor industry, material scientists, device physicists

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-Wave (PW) basis set
- LCAO basis set (Linear Combination of Atomic Orbitals) - *Note: separate module*
- NEGF (Non-Equilibrium Green's Functions)
- Pseudopotentials (Norm-conserving, PAW)
- Semi-empirical methods (DFTB, Slater-Koster)

## Capabilities
- **Plane-Wave DFT**:
  - Ground state energy and geometry optimization
  - Band structure and density of states
  - Stress and force calculations
  - Systematically convergent accuracy
- **Transport**:
  - I-V characteristics (via NEGF-LCAO/PW)
  - Transmission spectra
- **Properties**:
  - Optical properties
  - Phonons
  - Electron-phonon coupling
  - Mobility

## Key Strengths
- **Integrated Environment**: "NanoLab" GUI is one of the most advanced in the field.
- **Versatility**: Seamlessly switch between LCAO (fast/devices) and PlaneWave (accurate/bulk).
- **Industry Standard**: Widely used in semiconductor R&D (TCAD integration).
- **Python Scripting**: All functionality accessible via robust Python API.

## Computational Cost
- **Plane Wave**: Comparable to VASP/QE for similar cutoffs; $O(N^3)$.
- **LCAO**: $O(N)$ or $O(N^2)$ depending on method; faster for large systems.
- **License**: Commercial; cost is a significant factor.

## Best Practices
- **Basis Choice**: Use PlaneWave for bulk relaxation and high-precision electronic structure; use LCAO for large device/transport simulations.
- **Workflow**: Prototype with LCAO, verify with PW.
- **Scripting**: Use Python scripts for high-throughput workflows rather than just the GUI.

## Comparison with Other Codes
- **vs VASP/QE**: QuantumATK offers a much superior GUI and integrated transport (NEGF) workflows. VASP/QE are pure engines; QuantumATK is a platform.
- **vs SIESTA**: QuantumATK's LCAO mode is similar to SIESTA but integrated with a PW engine.

## Community and Support
- **Support**: Professional support from Synopsys.
- **Forum**: Active user forum.
- **Training**: Regular webinars and documented tutorials.

## Verification & Sources
**Primary sources**:
1. Synopsys QuantumATK Product Page: https://www.synopsys.com/silicon/quantumatk.html
2. QuantumATK Documentation: https://docs.quantumatk.com/
3. Smidstrup et al., J. Phys.: Condens. Matter 32, 015901 (2020).

**Confidence**: VERIFIED - Major commercial code.
https://www.synopsys.com/silicon/quantumatk.html
