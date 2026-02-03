# CALYPSO (Crystal structure AnaLYsis by Particle Swarm Optimization)

## Official Resources
- Homepage: http://www.calypso.cn/
- Documentation: http://www.calypso.cn/help/
- Source Repository: Distributed via website registration
- License: Proprietary (Free for academic use)

## Overview
CALYPSO is a software package for structure prediction using the Particle Swarm Optimization (PSO) algorithm. It can predict the stable structure of materials given only chemical composition and external conditions (pressure). CALYPSO is widely used for predicting crystal structures, 2D layers, clusters, and surfaces, and has been successful in predicting many high-pressure phases.

**Scientific domain**: Structure prediction, particle swarm optimization, materials discovery  
**Target user community**: Condensed matter physicists, materials scientists, chemists

## Theoretical Methods
- Particle Swarm Optimization (PSO)
- Structural evolution
- Symmetry constraints
- Local optimization via ab-initio codes
- Bond Characterization Matrix (BCM) for fingerprinting
- Artificial Bee Colony (ABC) algorithm (in newer versions)

## Capabilities (CRITICAL)
- Prediction of 3D crystals, 2D layers, 1D polymers, 0D clusters
- Surface structure prediction
- Interface structure prediction
- Variable composition search
- Superhard materials design
- Inverse design of materials with targeted properties (band gap, hardness)
- Integration with VASP, CASTEP, Quantum ESPRESSO, GULP, LAMMPS

**Sources**: CALYPSO website, Comp. Phys. Comm. 183, 1822 (2012)

## Inputs & Outputs
- **Input formats**: input.dat (parameters), split files for optimizer
- **Output data types**: Analysis results, predicted structures (POSCAR/CIF), evolutionary trajectory

## Interfaces & Ecosystem
- **Optimizers**: VASP, QE, CASTEP, etc.
- **Analysis**: CALYPSO_ANALYSIS tools
- **Visualization**: Compatible with standard structure viewers

## Workflow and Usage
1. Prepare `input.dat` with search parameters and system composition.
2. Prepare optimizer input files (e.g., `INCAR`, `POTCAR` for VASP).
3. Run CALYPSO executable.
4. CALYPSO generates structures -> Optimizer relaxes them -> CALYPSO updates swarm.
5. Analyze results using `calypso_analysis`.

## Performance Characteristics
- High search efficiency for small to medium systems
- Good global optimization capability
- Computationally expensive due to ab-initio relaxations

## Application Areas
- High-pressure physics
- Novel 2D materials
- Clusters and nano-objects
- Functional materials design

## Community and Support
- Large user community (thousands of licenses)
- Developed by Maosheng Miao and Yanming Ma groups (Jilin Univ./CSUN)
- Regular workshops
- Active forum

## Verification & Sources
**Primary sources**:
1. Homepage: http://www.calypso.cn/
2. Publication: Y. Wang, J. Lv, L. Zhu, and Y. Ma, Comp. Phys. Comm. 183, 1822 (2012)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: REGISTERED ACCESS
- Development: ACTIVE (Jilin University)
- Applications: PSO structure prediction, materials design, high pressure, 2D materials
