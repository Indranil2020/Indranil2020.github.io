# custodian

## Official Resources
- Homepage: https://materialsproject.github.io/custodian/
- Documentation: https://materialsproject.github.io/custodian/
- Source Repository: https://github.com/materialsproject/custodian
- License: MIT License

## Overview
Custodian is a simple, robust, and flexible just-in-time (JIT) job management framework in Python. It is designed to run scientific calculations (like VASP, Q-Chem, NwChem) while simultaneously monitoring them for errors and performing automatic error recovery. It acts as a "wrapper" around the executable, checking output files and logs in real-time to intervene if convergence fails or the job crashes.

**Scientific domain**: Error handling, job management, fault tolerance  
**Target user community**: Users of VASP, Q-Chem, and other stability-sensitive codes

## Capabilities (CRITICAL)
- **Error Detection**: Plugins for detecting specific errors in VASP (SCF divergence, walltime), Q-Chem, etc.
- **Error Recovery**: predefined "Handlers" to fix input files (e.g., changing mixing parameters, reducing time step) and restart the job.
- **Checkpointing**: Supports backup of files before runs.
- **Job Management**: Can run a sequence of jobs (e.g., relax -> static).
- **Extensible**: Easy to write new ErrorHandlers and Validators for any code.

**Sources**: Custodian documentation, Comp. Mater. Sci. 68, 314 (2013)

## Inputs & Outputs
- **Input formats**: Python script configuration
- **Output data types**: `custodian.json` report, corrected input files

## Interfaces & Ecosystem
- **VASP**: Extensive library of VASP error handlers
- **Q-Chem**: Support for Q-Chem errors
- **FireWorks/Atomate**: Used internally by these frameworks for reliability

## Workflow and Usage
1. Define jobs: `jobs = [VaspJob(output_file="vasp.out")]`
2. Define handlers: `handlers = [VaspErrorHandler(), UnconvergedErrorHandler()]`
3. Run custodian:
   ```python
   c = Custodian(handlers, jobs, max_errors=5)
   c.run()
   ```

## Performance Characteristics
- Lightweight wrapper
- Negligible overhead
- Significantly increases calculation success rate in high-throughput campaigns

## Application Areas
- High-throughput screening (running thousands of VASP jobs unsupervised)
- Long running simulations requiring stability

## Community and Support
- Developed by Materials Project team
- MIT License
- Widely used in the MP ecosystem

## Verification & Sources
**Primary sources**:
1. Homepage: https://materialsproject.github.io/custodian/
2. GitHub: https://github.com/materialsproject/custodian

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Error correction, VASP wrapper, reliability
