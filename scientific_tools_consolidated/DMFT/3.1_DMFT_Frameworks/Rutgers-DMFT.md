# Rutgers-DMFT

## Official Resources
- Homepage: https://github.com/HauleGroup/CODES
- Documentation: https://hauleweb.rutgers.edu/
- Source Repository: https://github.com/HauleGroup/CODES
- License: See repository for licensing details

## Overview
Rutgers-DMFT refers to the collection of DMFT codes developed by Kristjan Haule's group at Rutgers University, including EDMFTF and related impurity solvers. The repository contains various tools for DFT+DMFT calculations, continuous-time quantum Monte Carlo solvers, and analysis utilities developed and maintained by one of the leading DMFT research groups worldwide.

**Scientific domain**: DFT+DMFT, strongly correlated materials, quantum impurity solvers  
**Target user community**: Researchers performing advanced DFT+DMFT calculations on correlated materials

## Theoretical Methods
- DFT+DMFT (charge self-consistent)
- Continuous-time quantum Monte Carlo (CT-HYB, CT-QMC)
- Exact diagonalization
- Non-crossing approximation (NCA)
- Rotationally invariant master equation approach
- Anderson impurity model solvers
- Various computational methodologies from Haule group

## Capabilities (CRITICAL)
- Multiple impurity solver implementations
- DFT+DMFT workflows (EDMFTF)
- Advanced CTQMC solvers
- Analysis and post-processing tools
- Wien2k interface
- Wannier function tools
- Spectral function calculations
- Research-grade implementations
- Methods from cutting-edge DMFT research

**Sources**: Haule group repository (https://github.com/HauleGroup/CODES), Rutgers Physics, confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- DFT outputs (primarily Wien2k)
- Wannier projections
- Impurity model parameters
- DMFT configuration files

**Output data types**:
- Self-energies
- Green's functions
- Spectral functions
- Observables
- Convergence data

## Interfaces & Ecosystem
- **EDMFTF**: Primary DFT+DMFT framework
- **Wien2k**: Main DFT backend
- **Impurity solvers**: Multiple solver implementations
- **Rutgers Physics**: Leading DMFT research group
- **Haule group**: Kristjan Haule's research codes

## Limitations & Known Constraints
- Primarily research codes
- Documentation may be limited
- Requires expertise in DMFT methodology
- Installation and setup complex
- Some codes specific to Wien2k
- Platform: Linux/Unix
- Academic/research use focus

## Comparison with Other Codes
| Feature | Rutgers-DMFT (Haule) | TRIQS | w2dynamics |
| :--- | :--- | :--- | :--- |
| **Focus** | Materials / DFT+DMFT | General Framework / Solver Library | Multi-orbital Solvers |
| **DFT Interface** | Wien2k (Strongly integrated) | VASP, Wien2k, QE (via DFTTools) | VASP, Wien2k |
| **Solvers** | CT-QMC (Custom), NCA | CT-HYB, CT-INT, CT-SEG | CT-HYB, CT-INT |
| **Status** | Research / Academic Leader | Open Source Community Standard | Open Source |

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/HauleGroup/CODES
2. Haule group website: https://hauleweb.rutgers.edu/
3. Rutgers University Physics Department
4. EDMFTF documentation and papers

**Secondary sources**:
1. Haule group publications (>100 papers on DMFT)
2. EDMFTF users and tutorials
3. Rutgers strongly correlated materials research
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official repository: ACCESSIBLE (GitHub)
- Group website: ACCESSIBLE
- Source code: OPEN
- Developer: Kristjan Haule (leading DMFT expert)
- Research-grade codes
- Foundation for many DMFT developments
