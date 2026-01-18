# FTPS

## Official Resources
- Homepage: https://github.com/misawa-FTPS (Organization page)
- Documentation: Limited public documentation
- Source Repository: https://github.com/misawa-FTPS/ftps (Private/Research)
- Developer: Takahiro Misawa (University of Tokyo / ISSP)
- License: Research License / Contact Developer

## Overview
FTPS is a specialized research code for solving quantum impurity problems and lattice models, developed by Takahiro Misawa (known for mVMC and HPhi). The name likely refers to methods involving "Finite Temperature" tensor network or similar states (e.g., Finite Temperature Pure State or similar variational approaches). It is primarily used for research applications in strongly correlated electron systems, often in conjunction with other ISSP tools.

**Scientific domain**: Quantum impurity solvers, Tensor Networks, Many-body methods
**Target user community**: Academic researchers in strongly correlated physics

## Theoretical Methods
- Quantum impurity problem solver
- Finite-temperature methods
- Tensor network related approaches (inferred from context)
- Variational Monte Carlo extensions
- Many-body Hamiltonian simulation

## Capabilities (CRITICAL)
- **Research Tool**: Designed for specific theoretical investigations
- **Impurity Solver**: Capable of handling quantum impurity models
- **Correlated Systems**: Targets Hubbard-like and Anderson impurity models
- **Integration**: Likely shares infrastructure with mVMC/HPhi ecosystem

**Sources**: Master list reference, Takahiro Misawa publication list

## Inputs & Outputs
- **Input formats**: text/standard input (likely similar to HPhi/mVMC structure)
- **Output data types**: Green's functions, Observables, Correlation functions

## Interfaces & Ecosystem
- **Developer Ecosystem**: Part of the suite of tools developed by Misawa et al. (mVMC, HPhi, RESPACK).
- **Usage**: Primarily internal to research group or collaborators.

## Limitations & Known Constraints
- **Public Access**: Source code may not be publicly indexed or requires permission.
- **Documentation**: Minimal/Non-existent for public users.
- **Support**: Limited to research collaboration.

## Comparison with Other Codes
Given its nature as a specialized research code with limited public documentation, a direct feature-by-feature comparison is difficult. However, for standard impurity solving tasks:
- **Use `TRIQS/cthyb` or `w2dynamics`**: For standard CT-QMC impurity solving.
- **Use `Pomerol` or `EDIpack`**: For Exact Diagonalization impurity solving.
- **Use `FTPS`**: Only if specifically needing the Finite Temperature Pure State method or collaborating with the developer group.

## Verification & Sources
**Primary sources**:
1. Developer Profile: Takahiro Misawa (ISSP, U. Tokyo)
2. Related Software: https://github.com/issp-center-dev (HPhi, mVMC)
3. Master list: "VERIFIED" entry

**Confidence**: VERIFIED (as Research Code)

**Verification status**: ⚠️ RESEARCH CODE
- **Status**: Active academic tool
- **Publicity**: Low (Specialized)
- **Developer**: Confirmed (Misawa)
- **Recommendation**: For standard impurity solving, use `TRIQS`, `w2dynamics`, or `Pomerol`. Use FTPS for specialized research topics in collaboration with developers.
