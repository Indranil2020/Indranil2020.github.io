# AiiDA (Automated Interactive Infrastructure and Database for Computational Science)

## Official Resources
- Homepage: https://www.aiida.net/
- Documentation: https://aiida.readthedocs.io/
- Source Repository: https://github.com/aiidateam/aiida-core
- License: MIT License

## Overview
AiiDA is a sophisticated open-source workflow management system designed for high-throughput computational science. It focuses on reproducibility, provenance tracking, and automation. AiiDA manages the entire lifecycle of a calculation: preparing inputs, submitting jobs to supercomputers (via SLURM, PBS, etc.), parsing outputs, and storing the entire provenance graph (inputs, outputs, and the code that produced them) in a database.

**Scientific domain**: Workflow management, high-throughput computing, data provenance  
**Target user community**: Computational materials scientists, chemists, physicists

## Capabilities (CRITICAL)
- **Provenance Tracking**: Automatically records the full history of every calculation (inputs, codes, outputs) in a directed acyclic graph (DAG).
- **Workflow Engine**: Robust engine for running complex, interdependent calculations (WorkChains) with error handling.
- **Plugin System**: Massive ecosystem of plugins for DFT codes (VASP, QE, CP2K, etc.), schedulers, and transport.
- **HPC Interface**: Abstraction layer for interacting with remote clusters and job schedulers.
- **Database**: efficient storage and query of calculation data (PostgreSQL backend).
- **Querying**: Powerful QueryBuilder for searching the provenance graph.

**Sources**: AiiDA website, Comp. Mater. Sci. 111, 218 (2016)

## Inputs & Outputs
- **Input formats**: Python script defining `Process` or `WorkChain` inputs (Nodes).
- **Output data types**: Database entries (Nodes), files (retrieved from remote), provenance graph.

## Interfaces & Ecosystem
- **Plugins**: `aiida-quantumespresso`, `aiida-vasp`, `aiida-cp2k`, `aiida-siesta`, `aiida-wannier90`, etc.
- **Materials Cloud**: AiiDA archives can be published directly to Materials Cloud.
- **AiiDA Lab**: Web-based interface (Jupyter) for running AiiDA workflows.

## Workflow and Usage
1. **Setup**: Configure computer (remote cluster) and code (e.g., `pw.x`).
2. **Launch**: Submit a calculation or workflow via Python script or shell (`verdi run script.py`).
   ```python
   from aiida.engine import submit
   from aiida.plugins import CalculationFactory
   PwCalculation = CalculationFactory('quantumespresso.pw')
   submit(PwCalculation, **inputs)
   ```
3. **Monitor**: Check status with `verdi process list`.
4. **Analyze**: Query results with `QueryBuilder` or inspect nodes.

## Performance Characteristics
- High throughput: Can manage thousands of concurrent calculations.
- Daemon-based: Runs in the background, managing job submission and retrieval.
- Database performance depends on size and tuning (PostgreSQL).

## Application Areas
- High-throughput materials screening
- Database generation (e.g., Materials Cloud)
- Complex multi-step simulations (e.g., phonon dispersions, GW)
- Reproducible research studies

## Community and Support
- Developed by EPFL (Theos group) and wider community
- Large, active ecosystem
- Regular tutorials and hackathons
- Discourse forum for support

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.aiida.net/
2. GitHub: https://github.com/aiidateam/aiida-core
3. Publication: G. Pizzi et al., Comp. Mater. Sci. 111, 218 (2016); S. P. Huber et al., Sci. Data 7, 300 (2020)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (EPFL/AiiDA Team)
- Applications: Workflow management, provenance, high-throughput
