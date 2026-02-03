# String methods

## Official Resources
- Homepage: Method - Implemented in many codes (VASP, Quantum ESPRESSO, etc.)
- Documentation: https://theory.cm.utexas.edu/vtsttools/ (VTST implementation)
- Source Repository: https://github.com/henkelmanlab/vtstscripts
- License: Varies by implementation

## Overview
String methods (including the Simplified String Method and Finite Temperature String Method) are a class of chain-of-states methods used to find the minimum energy path (MEP) or transition pathways in complex energy landscapes. Similar to NEB, they evolve a string of images connecting two minima, but with different parametrization and evolution equations, often offering better stability or different convergence properties.

**Scientific domain**: Transition path sampling, rare events, minimum energy paths  
**Target user community**: Computational chemists, condensed matter physicists

## Theoretical Methods
- String Method (SM)
- Simplified String Method (SSM)
- Finite Temperature String Method (FTSM)
- Zero-Temperature String Method
- Reparametrization of path
- Transition Path Theory

## Capabilities (CRITICAL)
- Finding Minimum Energy Paths (MEP)
- Determining transition tubes in free energy landscapes
- Sampling rare events
- Calculating transition rates
- Handling collective variables
- Implemented in: VASP (VTST), Quantum ESPRESSO, molecular dynamics codes

**Sources**: E, Ren, and Vanden-Eijnden, Phys. Rev. B 66, 052301 (2002)

## Inputs & Outputs
- **Input formats**: Initial/Final states, string of images
- **Output data types**: Evolved path, energy profile, free energy surface

## Interfaces & Ecosystem
- **VASP**: Via VTST tools
- **Quantum ESPRESSO**: Native or plugin support
- **MD Codes**: Often implemented for exploring free energy surfaces

## Workflow and Usage
1. Define collective variables or reaction coordinate
2. Initialize string connecting reactants and products
3. Evolve string according to potential/free energy gradient
4. Reparametrize string to maintain equidistant images
5. Converge to MEP or principal curve

## Performance Characteristics
- Highly parallelizable (images evolved independently step-wise)
- Convergence depends on landscape complexity
- FTSM computationally expensive due to sampling

## Application Areas
- Phase transitions
- Protein conformational changes
- Nucleation processes
- Chemical reaction pathways
- Diffusion in complex media

## Community and Support
- Method developers (E, Vanden-Eijnden)
- Code-specific communities (VASP, QE)
- Literature and theoretical physics community

## Verification & Sources
**Primary sources**:
1. Publication: W. E, W. Ren, and E. Vanden-Eijnden, Phys. Rev. B 66, 052301 (2002)
2. VTST Tools: https://theory.cm.utexas.edu/vtsttools/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Method: STANDARD
- Documentation: Available in implementation codes
- Applications: MEP finding, rare events, transition paths
