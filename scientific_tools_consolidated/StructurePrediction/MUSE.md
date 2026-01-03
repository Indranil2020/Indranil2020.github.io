# MUSE (Multi-algorithm collaborative Universal Structure-prediction Environment)

## Official Resources
- **Homepage**: https://github.com/zhongliliu/muse
- **Source Repository**: https://github.com/zhongliliu/muse
- **Documentation**: https://github.com/zhongliliu/muse/tree/master/doc
- **License**: GPL-3.0

## Overview
MUSE (Multi-algorithm collaborative Universal Structure-prediction Environment) is a crystal structure prediction package that combines multiple global optimization algorithms—Evolutionary Algorithm, Simulated Annealing, and Basin Hopping—to efficiently locate stable structures. It is designed to leverage the strengths of different search strategies.

**Scientific domain**: Crystal structure prediction, global optimization  
**Target user community**: Materials physicists, condensed matter researchers

## Theoretical Methods
- **Hybrid Optimization**: Collaborative use of:
  - Evolutionary Algorithm (EA)
  - Simulated Annealing (SA)
  - Basin Hopping (BH)
- **Density Functional Theory**: Interfaces with DFT codes for energy evaluation.
- **Symmetry Constraints**: Uses symmetry to reduce search space.

## Capabilities
- **Crystal Structure Prediction**: For 3D crystals under pressure.
- **Collaborative Search**: different algorithms run in parallel and exchange information.
- **Interfaces**: Supports VASP, CASTEP, Quantum ESPRESSO, GULP.
- **Pressure Support**: High-pressure phase prediction.

## Inputs & Outputs
- **Input formats**: `input.dat` (control parameters), `seeds` directory.
- **Output data types**: `low_energy_structures`, optimization logs.

## Interfaces & Ecosystem
- **VASP**: Primary DFT interface.
- **CASTEP/QE**: Supported interfaces.
- **GULP**: For empirical potentials.

## Verification & Sources
- **Confidence**: ✅ VERIFIED
- **Primary Source**: [MUSE GitHub](https://github.com/zhongliliu/muse)
- **Reference**: Z. L. Liu, "MUSE: Multi-algorithm collaborative crystal structure prediction", *Comput. Phys. Commun.* 185, 1893 (2014).
