# Materials Cloud

## Official Resources
- Homepage: https://www.materialscloud.org/
- Documentation: https://www.materialscloud.org/learn
- Source Repository: https://github.com/materialscloud-org
- License: Various (platform components are open source)

## Overview
Materials Cloud is a web platform for Open Science in computational materials science. Built on top of AiiDA, it enables researchers to FAIRly share their data, calculations, and provenance graphs. It provides interactive tools for visualizing data (structures, bands, phonons) and a "Work" section for running simulations in the cloud (AiiDA Lab).

**Scientific domain**: Open science, data repository, interactive visualization  
**Target user community**: Computational materials scientists, AiiDA users

## Capabilities (CRITICAL)
- **Archive**: Long-term repository for research data with DOI assignment.
- **Discover**: Interactive exploration of curated datasets (e.g., 2D materials, COFs).
- **Explore**: Graphical interface to browse AiiDA provenance graphs.
- **Tools**: Web-based tools for common tasks (seekpath, band structure plotting, structure visualization).
- **AiiDA Lab**: Cloud environment (JupyterHub) for running AiiDA workflows.

**Sources**: Materials Cloud website, APL Mater. 6, 101101 (2018)

## Inputs & Outputs
- **Input formats**: AiiDA export files, various structure formats
- **Output data types**: DOIs, interactive plots, downloadable datasets

## Interfaces & Ecosystem
- **AiiDA**: Deeply integrated; acts as the frontend for AiiDA data.
- **Quantum Mobile**: Virtual machine with pre-installed codes.
- **OPTIMADE**: Supports the OPTIMADE API for interoperability.

## Workflow and Usage
1. Run calculations with AiiDA.
2. Export database: `verdi archive create data.aiida`.
3. Upload to Materials Cloud Archive.
4. Data becomes browsable and citable.

## Performance Characteristics
- Web-based; performance depends on browser and server load.
- Visualizers are highly optimized for interactivity.

## Application Areas
- Publishing reproducible research data
- Teaching (via AiiDA Lab)
- Quick visualization of standard file formats

## Community and Support
- Developed by EPFL (Theos), PSI, and partners
- Core component of the MARVEL NCCR

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.materialscloud.org/
2. Publication: N. Mounet et al., Nat. Nanotechnol. 13, 246 (2018); A. V. Yakutovich et al., Comp. Mater. Sci. 188, 110165 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Data sharing, visualization, AiiDA integration
