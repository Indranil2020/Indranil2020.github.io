# ASR (Atomic Simulation Recipes)

## Official Resources
- Homepage: https://asr.readthedocs.io/
- Documentation: https://asr.readthedocs.io/
- Source Repository: https://gitlab.com/ase/asr
- License: GNU General Public License v3.0

## Overview
Atomic Simulation Recipes (ASR) is a Python framework for defining and executing simulation workflows. Developed at DTU, it allows users to define "recipes" (workflows) that link calculations (typically using ASE and GPAW) to results. ASR is the engine behind the C2DB database. It emphasizes caching, reproducibility, and automatic presentation of results via a web interface.

**Scientific domain**: Workflow automation, database generation, reproducibility  
**Target user community**: ASE/GPAW users, C2DB contributors

## Capabilities (CRITICAL)
- **Recipes**: Modular Python scripts defining calculations (e.g., `relax.py`, `bandstructure.py`).
- **Caching**: Smart caching of results ("smart make" functionality) to avoid recomputing.
- **Web UI**: Automatically generates a web page for browsing results in a directory.
- **CLI**: Command line interface (`asr run relax`).
- **Dependency**: Automatic handling of task dependencies.

**Sources**: ASR documentation, Adv. Theory Simul. 2, 1900080 (2019)

## Inputs & Outputs
- **Input formats**: ASE atoms objects, CLI arguments
- **Output data types**: Pickled result files (`results-relax.pkl`), JSON

## Interfaces & Ecosystem
- **ASE**: Core object model.
- **GPAW**: Primary calculator supported (but code agnostic in principle).
- **MyQueue**: Scheduler integration.

## Workflow and Usage
1. Create a directory for a material.
2. Run recipes:
   ```bash
   asr run relax
   asr run bandstructure
   ```
3. View results: `asr browser` starts a local web server.

## Performance Characteristics
- Very low overhead
- Excellent for managing "high-throughput" studies of a moderate number of materials (100s-1000s)

## Application Areas
- Database generation (C2DB)
- Systematic material studies
- Reproducible thesis work

## Community and Support
- Developed by CAMD (DTU)
- Active development

## Verification & Sources
**Primary sources**:
1. Documentation: https://asr.readthedocs.io/
2. GitLab: https://gitlab.com/ase/asr

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitLab)
- Development: ACTIVE
- Applications: Workflows, recipes, caching
