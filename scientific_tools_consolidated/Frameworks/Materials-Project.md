# Materials Project

## Official Resources
- Homepage: https://materialsproject.org/
- Documentation: https://docs.materialsproject.org/
- Source Repository: https://github.com/materialsproject
- License: CC BY 4.0 (Data), MIT (Code)

## Overview
The Materials Project (MP) is a major initiative to compute the properties of all known inorganic materials and provide the data freely to the public. It provides a searchable database of calculated properties (band structures, elastic constants, piezoelectricity, etc.) generated using high-throughput DFT (VASP). It also provides an API (MAPIDoc) and a Python client (`mp_api`) for programmatic access.

**Scientific domain**: Materials database, high-throughput DFT, materials informatics  
**Target user community**: Materials scientists, data scientists, chemists

## Capabilities (CRITICAL)
- **Data Access**: Web interface to search materials by elements, formula, or properties.
- **Properties**: Structure, formation energy, band gap, density of states, elastic tensor, piezoelectric tensor, dielectric constant, X-ray absorption spectra, phonon dispersion.
- **Apps**: Phase diagram app, Pourbaix diagram app, Battery explorer, Crystal toolkit.
- **API**: RESTful API (`MPRester`) for downloading large datasets.
- **Provenance**: Full calculation details (inputs/outputs) available.

**Sources**: Materials Project website, APL Mater. 1, 011002 (2013)

## Inputs & Outputs
- **Input formats**: Queries (JSON, Python criteria)
- **Output data types**: Structures (CIF/POSCAR), Property values, Plots

## Interfaces & Ecosystem
- **Pymatgen**: The official analysis library and API client.
- **FireWorks/Atomate**: The workflow stack used to generate the data.
- **Crystal Toolkit**: Web-based visualization.

## Workflow and Usage
1. **Web**: Go to materialsproject.org, search for "Li-Fe-P", filter by stable materials.
2. **Python**:
   ```python
   from mp_api.client import MPRester
   with MPRester("API_KEY") as mpr:
       docs = mpr.summary.search(chemsys="Li-Fe-P", band_gap=(1.0, 3.0))
   ```

## Performance Characteristics
- Massive database (>140,000 materials)
- API rate limits apply
- Pre-computed data (instant retrieval vs waiting for calculation)

## Application Areas
- Battery material discovery
- Photovoltaic screening
- Thermoelectric material search
- Machine learning model training

## Community and Support
- Developed by LBNL and partners
- Huge user base (>200,000 registered users)
- Active forum (MatSci.org)

## Verification & Sources
**Primary sources**:
1. Homepage: https://materialsproject.org/
2. Documentation: https://docs.materialsproject.org/
3. Publication: A. Jain et al., APL Mater. 1, 011002 (2013)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (Code), OPEN ACCESS (Data)
- Development: ACTIVE
- Applications: Materials database, API, discovery
