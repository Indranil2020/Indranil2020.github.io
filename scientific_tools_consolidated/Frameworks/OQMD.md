# OQMD (Open Quantum Materials Database)

## Official Resources
- Homepage: https://oqmd.org/
- Documentation: https://oqmd.org/static/docs/index.html
- Source Repository: https://github.com/wolverton-research-group/qmpy (Software stack)
- License: CC BY 4.0 (Data)

## Overview
The Open Quantum Materials Database (OQMD) is a high-throughput database of DFT-calculated thermodynamic and structural properties of materials. It focuses heavily on thermodynamics, phase stability, and the discovery of new stable compounds. It is built using the `qmpy` software stack.

**Scientific domain**: Materials database, thermodynamics, phase stability  
**Target user community**: Materials scientists, metallurgists

## Capabilities (CRITICAL)
- **Database**: Contains >1 million materials (calculated + experimental).
- **Phase Diagrams**: Grand canonical linear programming for constructing convex hulls and phase diagrams.
- **Stability**: Formation energies and stability analysis (energy above hull).
- **Properties**: Crystal structure, formation energy, band gap, magnetic moment.
- **Software**: `qmpy` (Quantum Materials Python) handles the workflow and analysis.

**Sources**: OQMD website, JOM 65, 1501 (2013)

## Inputs & Outputs
- **Input formats**: Query via web or API
- **Output data types**: Structure files, phase diagrams, thermodynamic data

## Interfaces & Ecosystem
- **VASP**: Calculation engine used for the database
- **qmpy**: The Python framework powering OQMD
- **API**: REST API for data access

## Workflow and Usage
1. **Web**: Search for "Fe-O" to see the phase diagram and stable compounds.
2. **API**: Use `qmpy` or standard HTTP requests to download formation energies.
   ```python
   import qmpy
   entries = qmpy.Entry.objects.filter(element_set="Fe-O")
   ```

## Performance Characteristics
- Large database focused on ground state thermodynamics
- `qmpy` automates the VASP workflows

## Application Areas
- Discovery of new stable phases
- Battery materials (voltage profiles)
- Alloy design

## Community and Support
- Developed by Wolverton Group (Northwestern University)
- Open data policy

## Verification & Sources
**Primary sources**:
1. Homepage: https://oqmd.org/
2. Publication: J. E. Saal et al., JOM 65, 1501 (2013)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (qmpy on GitHub)
- Development: ACTIVE
- Applications: Materials database, thermodynamics
