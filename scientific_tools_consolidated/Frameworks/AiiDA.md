# AiiDA

## Official Resources
- Homepage: https://www.aiida.net/
- Documentation: https://aiida.readthedocs.io/
- Source Repository: https://github.com/aiidateam/aiida-core
- License: MIT License

## Overview
AiiDA (Automated Interactive Infrastructure and Database for computational science) is a Python framework for computational science workflow management and provenance tracking. Developed at EPFL, it automatically preserves the full provenance of all calculations, enabling reproducibility, data sharing, and advanced queries. AiiDA is designed for high-throughput materials science but applicable to any computational workflow requiring data management and automation.

**Scientific domain**: Workflow management, provenance tracking, computational infrastructure  
**Target user community**: Computational scientists, high-throughput researchers, database developers

## Theoretical Methods
AiiDA itself does not implement methods but provides infrastructure for:
- Workflow automation
- Provenance tracking (directed acyclic graph)
- Database management (PostgreSQL)
- Plugin system for codes
- Remote computation management
- Data querying and analysis

## Capabilities (CRITICAL)
- Automated workflow management and execution
- Full data provenance tracking (inputs, outputs, code versions)
- Reproducible calculations
- High-throughput materials screening
- Job submission to local and remote clusters
- RESTful API for data access
- Powerful query language for database
- Workflow error handling and recovery
- Plugin system for any computational code
- Caching to avoid redundant calculations
- Share and export calculations (AiiDA archives)
- Integration with Materials Cloud
- Python API for workflow development
- Graph-based workflow representation

**Sources**: Official AiiDA documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Python scripts defining workflows
  - Calculation inputs managed through plugins
  - Structure files (CIF, XYZ, ASE, pymatgen)
  
- **Output data types**:
  - PostgreSQL database storing all data
  - AiiDA archive files (.aiida)
  - JSON export format
  - Calculation outputs accessed via Python API
  - Provenance graphs

## Interfaces & Ecosystem
- **Plugin ecosystem**:
  - aiida-quantumespresso (Quantum ESPRESSO)
  - aiida-vasp (VASP)
  - aiida-cp2k (CP2K)
  - aiida-fleur (Fleur)
  - aiida-wannier90 (Wannier90)
  - aiida-yambo (Yambo)
  - 50+ official and community plugins
  
- **Materials databases**:
  - Materials Cloud integration
  - OPTIMADE API support
  - Export to external databases
  
- **Workflow engines**:
  - Built-in workflow engine
  - Process management and scheduling
  - Daemon for background execution

## Limitations & Known Constraints
- **Learning curve**: Steep; requires Python and workflow concepts
- **Database setup**: PostgreSQL and RabbitMQ required
- **Plugin dependency**: Calculation codes need plugins (most major codes covered)
- **Computational overhead**: Framework overhead for small calculations
- **Migration**: Database schema changes between versions
- **Documentation**: Extensive but can be overwhelming
- **Platform**: Primarily Linux/Unix; macOS supported

## Verification & Sources
**Primary sources**:
1. Official website: https://www.aiida.net/
2. Documentation: https://aiida.readthedocs.io/
3. GitHub repository: https://github.com/aiidateam/aiida-core
4. S. P. Huber et al., Sci. Data 7, 300 (2020) - AiiDA 1.0
5. G. Pizzi et al., Comput. Mater. Sci. 111, 218 (2016) - AiiDA original paper

**Secondary sources**:
1. AiiDA tutorials and workshops
2. Materials Cloud documentation
3. Plugin registry
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Very active (Slack, forum, workshops)
- Academic citations: >300
- Active development: Continuous releases, MARVEL funded
