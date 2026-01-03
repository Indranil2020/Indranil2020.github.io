# JARVIS (Joint Automated Repository for Various Integrated Simulations)

## Official Resources
- Homepage: https://jarvis.nist.gov/
- Documentation: https://jarvis-tools.readthedocs.io/
- Source Repository: https://github.com/usnistgov/jarvis
- License: NIST License (Public Domain)

## Overview
JARVIS is an integrated framework and database developed by NIST for data-driven materials design. It consists of the JARVIS-DFT database (DFT calculations), JARVIS-FF (Force fields), and JARVIS-ML (Machine learning). It emphasizes properties relevant to applications (e.g., solar cells, thermoelectrics, dielectrics) and uses high-throughput workflows powered by `jarvis-tools`.

**Scientific domain**: Materials database, high-throughput DFT, machine learning  
**Target user community**: Materials scientists, ML researchers

## Capabilities (CRITICAL)
- **JARVIS-DFT**: Database of >80,000 materials with VASP/TB-mBJ calculations.
- **Properties**: Elastic constants, dielectric tensors, piezoelectricity, topological invariants, solar efficiency, exfoliation energy.
- **JARVIS-FF**: Database of classical force field properties.
- **JARVIS-ML**: Machine learning models and descriptors.
- **Tools**: `jarvis-tools` python package for automation and analysis.

**Sources**: JARVIS website, Sci. Data 5, 180082 (2018)

## Inputs & Outputs
- **Input formats**: VASP inputs, LAMMPS inputs
- **Output data types**: JSON databases, XML files

## Interfaces & Ecosystem
- **VASP**: Primary DFT engine
- **LAMMPS**: Primary MD engine
- **TensorFlow/PyTorch**: Used for ML models
- **Web API**: For programmatic access

## Workflow and Usage
1. **Web**: Search JARVIS-DFT for "MoS2".
2. **Python**:
   ```python
   from jarvis.db.figshare import data
   dft_data = data('dft_3d')
   ```

## Performance Characteristics
- High-quality data (TB-mBJ for band gaps)
- Comprehensive coverage of 2D materials

## Application Areas
- 2D materials discovery
- Thermoelectric materials
- Topological materials
- Machine learning benchmarking

## Community and Support
- Developed by NIST (Kamal Choudhary et al.)
- Public domain

## Verification & Sources
**Primary sources**:
1. Homepage: https://jarvis.nist.gov/
2. GitHub: https://github.com/usnistgov/jarvis
3. Publication: K. Choudhary et al., Sci. Data 5, 180082 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (NIST)
- Applications: Materials database, ML, DFT
