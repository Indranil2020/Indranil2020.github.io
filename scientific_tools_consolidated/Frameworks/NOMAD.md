# NOMAD (Novel Materials Discovery)

## Official Resources
- Homepage: https://nomad-lab.eu/
- Documentation: https://nomad-lab.eu/prod/v1/staging/docs/
- Source Repository: https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR
- License: Apache License 2.0

## Overview
NOMAD is a web-based data management platform for materials science. Unlike MP or OQMD which generate their own data, NOMAD is primarily a repository for *archiving* and *sharing* data produced by the community from *any* code. It processes uploaded raw output files (from VASP, QE, etc.), normalizes them into a code-independent format (NOMAD Metainfo), and makes them searchable and analyzable.

**Scientific domain**: Materials data repository, FAIR data, materials informatics  
**Target user community**: All computational materials scientists

## Capabilities (CRITICAL)
- **Repository**: Long-term archival of raw calculation data (upload inputs/outputs).
- **Normalization**: Parsers for >40 electronic structure codes convert data to a unified schema.
- **Encyclopedia**: User-friendly view of materials properties derived from the archive.
- **AI Toolkit**: Jupyter-based environment for analyzing NOMAD data using machine learning.
- **ELN**: Electronic Lab Notebook features for managing experimental and computational workflows.

**Sources**: NOMAD website, Nature 604, 635 (2022)

## Inputs & Outputs
- **Input formats**: Raw output files from almost any electronic structure code (VASP, QE, Abinit, Crystal, etc.)
- **Output data types**: Archived datasets, normalized JSON/metainfo

## Interfaces & Ecosystem
- **Supported Codes**: VASP, Quantum ESPRESSO, Abinit, FHI-aims, Exciting, CP2K, Gaussian, etc.
- **API**: REST API for querying the archive.
- **Python Library**: `nomad-lab` package for interacting with the API.

## Workflow and Usage
1. Run calculations locally.
2. Zip the directory.
3. Upload to NOMAD Repository (publish or keep private).
4. NOMAD automatically parses and indexes the data.
5. Others can find and download your raw files for reproducibility.

## Performance Characteristics
- Largest database of materials calculations in the world (>100 million entries)
- Focus on raw data (not just derived properties)

## Application Areas
- Reproducibility (access to exact inputs/outputs)
- Big data analytics / Machine Learning
- Benchmarking codes

## Community and Support
- European Center of Excellence (Fairmat)
- Large international consortium
- Hosted by MPCDF

## Verification & Sources
**Primary sources**:
1. Homepage: https://nomad-lab.eu/
2. Publication: C. Draxl and M. Scheffler, MRS Bulletin 43, 676 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitLab)
- Development: ACTIVE
- Applications: Data repository, FAIR data, parsing
