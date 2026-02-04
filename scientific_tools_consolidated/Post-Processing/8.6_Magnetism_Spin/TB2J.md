# TB2J

## Official Resources
- Homepage: https://tb2j.readthedocs.io/
- Documentation: https://tb2j.readthedocs.io/en/latest/
- Source Repository: https://github.com/mailhexu/TB2J
- License: BSD 2-Clause License

## Overview
TB2J is a Python package for calculating magnetic interaction parameters (isotropic exchange J, Dzyaloshinskii-Moriya interaction D, and magnetic anisotropy K) from first-principles calculations. It uses the magnetic force theorem and Green's function method within a tight-binding representation (Wannier functions or LCAO). TB2J interfaces with various DFT codes and outputs parameters for spin dynamics codes like Spirit, VAMPIRE, and MultiBader.

**Scientific domain**: Magnetism, exchange interactions, spintronics  
**Target user community**: Magnetism researchers, DFT users

## Theoretical Methods
- Magnetic Force Theorem (Liechtenstein formula)
- Green's function method
- Local Force Theorem
- Wannier function projection (Wannier90)
- LCAO basis (SIESTA, OpenMX)
- Heisenberg Hamiltonian mapping

## Capabilities (CRITICAL)
- Calculation of isotropic exchange (J_ij)
- Calculation of Dzyaloshinskii-Moriya interaction (DMI)
- Calculation of single-ion anisotropy (SIA)
- Orbital-resolved exchange interactions
- Support for collinear and non-collinear magnetism
- Spin-orbit coupling effects
- Automatic generation of input for spin dynamics codes

**Sources**: TB2J documentation, Comp. Phys. Comm. 261, 107760 (2021)

## Inputs & Outputs
- **Input formats**: Wannier90 output (wannier90_hr.dat), SIESTA/OpenMX Hamiltonian files
- **Output data types**: Exchange parameters (txt, xml), geometrical structure, plotting scripts

## Interfaces & Ecosystem
- **DFT Codes**: Wannier90 (VASP, QE, etc.), SIESTA, OpenMX, ABINIT
- **Spin Dynamics**: Generates inputs for Spirit, VAMPIRE, UppASD
- **Multibinit**: Interface available

## Workflow and Usage
1. Perform DFT calculation (ground state).
2. Generate Wannier functions (Wannier90) or use LCAO output (SIESTA).
3. Run TB2J: `wannier2tb2j.py` or `siesta2tb2j.py`.
4. Analyze output `exchange.out` or use generated files for spin dynamics.

## Performance Characteristics
- Very efficient post-processing (Green's function integration)
- Memory usage depends on tight-binding basis size
- Parallelization over k-points

## Application Areas
- Magnon dispersion relations
- Curie temperature estimation
- Skyrmion stabilization mechanisms (DMI)
- 2D magnetic materials
- Complex magnetic textures

## Community and Support
- Open-source (BSD)
- Developed by He Xu (Liège/Luxembourg)
- Active GitHub repository and documentation

## Verification & Sources
**Primary sources**:
1. Homepage: https://tb2j.readthedocs.io/
2. GitHub: https://github.com/mailhexu/TB2J
3. Publication: X. He et al., Comp. Phys. Comm. 266, 108036 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Magnetic exchange, DMI, DFT interface, spin dynamics
