#!/usr/bin/env python3
"""
Generate ALL remaining skeleton .md files for 372 tools
Complete coverage of all categories
"""

import os
from pathlib import Path

# Complete tool list for remaining categories
remaining_tools = {
    "Tight-Binding": [
        ("WannierTools", "CONFIRMED", "https://github.com/quanshengwu/wannier_tools"),
        ("WannierBerri", "VERIFIED", "UNKNOWN"),
        ("pythtb", "CONFIRMED", "http://www.physics.rutgers.edu/pythtb/"),
        ("TBmodels", "CONFIRMED", "https://tbmodels.greschd.ch/"),
        ("Z2Pack", "CONFIRMED", "UNKNOWN"),
        ("Kwant", "VERIFIED", "UNKNOWN"),
        ("Pybinding", "LOW_CONF", "UNKNOWN"),
        ("TBSTUDIO", "LOW_CONF", "UNKNOWN"),
        ("TopoTB", "LOW_CONF", "https://github.com/xlhuang-phy/TopoTB"),
        ("TBPLaS", "LOW_CONF", "https://github.com/tbplas/tbplas"),
        ("Chinook", "LOW_CONF", "UNKNOWN"),
        ("BoltzWann", "VERIFIED", "UNKNOWN"),
        ("PyWannier90", "LOW_CONF", "UNKNOWN"),
        ("WOPT", "LOW_CONF", "UNKNOWN"),
        ("VASP2Wannier90", "LOW_CONF", "UNKNOWN"),
        ("ir2tb", "VERIFIED", "https://github.com/yuzie007/ir2tb"),
        ("RESPACK", "LOW_CONF", "https://respack.org/"),
        ("TightBinding++", "LOW_CONF", "UNKNOWN"),
        ("QuantumLattice", "UNCERTAIN", "UNKNOWN"),
        ("QuantNBody", "LOW_CONF", "https://quantnbody.readthedocs.io/"),
        ("Paoflow", "LOW_CONF", "UNKNOWN"),
        ("MagneticTB", "LOW_CONF", "UNKNOWN"),
        ("MagneticKP", "LOW_CONF", "UNKNOWN"),
        ("FermiSurfer", "CONFIRMED", "https://fermisurfer.osdn.jp/"),
    ],
    
    "Phonons": [
        ("ShengBTE", "CONFIRMED", "http://www.shengbte.org/"),
        ("ALAMODE", "CONFIRMED", "https://alamode.readthedocs.io/"),
        ("almaBTE", "CONFIRMED", "https://www.almabte.org/"),
        ("TDEP", "CONFIRMED", "https://ollehellman.github.io/"),
        ("EPW", "CONFIRMED", "https://docs.epw-code.org/"),
        ("PERTURBO", "VERIFIED", "UNKNOWN"),
        ("Phoebe", "VERIFIED", "UNKNOWN"),
        ("PHON", "LOW_CONF", "UNKNOWN"),
        ("PHONON", "LOW_CONF", "UNKNOWN"),
        ("YPHON", "LOW_CONF", "UNKNOWN"),
        ("ATAT", "LOW_CONF", "UNKNOWN"),
        ("FROPHO", "LOW_CONF", "UNKNOWN"),
        ("hiPhive", "LOW_CONF", "UNKNOWN"),
        ("ASE-phonons", "UNCERTAIN", "UNKNOWN"),
        ("kALDo", "LOW_CONF", "UNKNOWN"),
        ("GPU_PBTE", "LOW_CONF", "UNKNOWN"),
        ("PhonTS", "UNCERTAIN", "UNKNOWN"),
        ("SCAILD", "LOW_CONF", "UNKNOWN"),
        ("QSCAILD", "LOW_CONF", "UNKNOWN"),
        ("SSCHA", "LOW_CONF", "UNKNOWN"),
        ("ALM", "LOW_CONF", "UNKNOWN"),
        ("thirdorder.py", "LOW_CONF", "UNKNOWN"),
        ("THERMACOND", "LOW_CONF", "UNKNOWN"),
        ("ALATDYN", "LOW_CONF", "UNKNOWN"),
        ("OpenBTE", "LOW_CONF", "UNKNOWN"),
        ("DMDW", "LOW_CONF", "UNKNOWN"),
        ("RTDW", "LOW_CONF", "UNKNOWN"),
        ("epiq", "LOW_CONF", "UNKNOWN"),
        ("API_Phonons", "LOW_CONF", "UNKNOWN"),
        ("Phonopy-API", "LOW_CONF", "UNKNOWN"),
        ("Pheasy", "LOW_CONF", "UNKNOWN"),
        ("Simphony", "LOW_CONF", "UNKNOWN"),
    ],
    
    "Dynamics": [
        ("i-PI", "CONFIRMED", "https://ipi-code.org/"),
        ("LAMMPS", "VERIFIED", "https://www.lammps.org/"),
        ("PLUMED", "VERIFIED", "UNKNOWN"),
        ("GROMACS", "VERIFIED", "UNKNOWN"),
        ("AMBER", "LOW_CONF", "UNKNOWN"),
        ("CHARMM", "LOW_CONF", "UNKNOWN"),
        ("NAMD", "LOW_CONF", "UNKNOWN"),
        ("DL_POLY", "LOW_CONF", "UNKNOWN"),
        ("N2P2", "LOW_CONF", "UNKNOWN"),
        ("DeepMD-kit", "LOW_CONF", "UNKNOWN"),
        ("OpenMD", "LOW_CONF", "https://openmd.org/"),
        ("IMD", "LOW_CONF", "UNKNOWN"),
        ("NEB", "VERIFIED", "UNKNOWN"),
        ("String-methods", "LOW_CONF", "UNKNOWN"),
        ("Metadynamics", "LOW_CONF", "UNKNOWN"),
        ("libAtoms-Quippy", "LOW_CONF", "UNKNOWN"),
        ("MDI-MolSSI", "LOW_CONF", "UNKNOWN"),
    ],
    
    "Structure-Prediction": [
        ("XtalOpt", "CONFIRMED", "https://xtalopt.github.io/"),
        ("CALYPSO", "CONFIRMED", "http://www.calypso.cn/"),
        ("AIRSS", "CONFIRMED", "https://www.mtg.msm.cam.ac.uk/Codes/AIRSS"),
        ("GASP", "CONFIRMED", "https://github.com/iwan-w/gasp"),
        ("MAISE", "LOW_CONF", "UNKNOWN"),
        ("EVO", "LOW_CONF", "UNKNOWN"),
        ("FLAME", "LOW_CONF", "UNKNOWN"),
        ("Basin-hopping", "LOW_CONF", "UNKNOWN"),
        ("HTOCSP", "LOW_CONF", "UNKNOWN"),
        ("PyXtal", "LOW_CONF", "UNKNOWN"),
        ("PXRDGen", "LOW_CONF", "UNKNOWN"),
        ("OpenCSP", "LOW_CONF", "UNKNOWN"),
        ("GMIN", "LOW_CONF", "UNKNOWN"),
        ("ASE-GA", "UNCERTAIN", "UNKNOWN"),
        ("ASE-BasinHopping", "LOW_CONF", "UNKNOWN"),
        ("MUSE", "UNCERTAIN", "UNKNOWN"),
        ("PyMaterial-Search", "UNCERTAIN", "UNKNOWN"),
        ("PyMetadynamics", "UNCERTAIN", "UNKNOWN"),
        ("MaterialsProject-ML", "LOW_CONF", "UNKNOWN"),
        ("PyXtal-ML", "LOW_CONF", "UNKNOWN"),
        ("Oganov-ML", "UNCERTAIN", "UNKNOWN"),
    ],
    
    "Post-Processing": [
        ("vaspkit", "CONFIRMED", "https://vaspkit.com/"),
        ("sumo", "CONFIRMED", "https://sumo.readthedocs.io/"),
        ("pyprocar", "CONFIRMED", "https://pyprocar.readthedocs.io/"),
        ("PyARPES", "LOW_CONF", "UNKNOWN"),
        ("BandUP", "LOW_CONF", "https://github.com/bandup/bandup"),
        ("fold2Bloch", "LOW_CONF", "UNKNOWN"),
        ("irvsp", "LOW_CONF", "UNKNOWN"),
        ("SeeK-path", "LOW_CONF", "UNKNOWN"),
        ("PyProcar-Unfold", "UNCERTAIN", "UNKNOWN"),
        ("IrRep", "LOW_CONF", "UNKNOWN"),
        ("effectivemass", "LOW_CONF", "UNKNOWN"),
        ("BerryPI", "LOW_CONF", "UNKNOWN"),
        ("Chern-Number", "LOW_CONF", "UNKNOWN"),
        ("Berry-Phase", "LOW_CONF", "UNKNOWN"),
        ("BoltzTraP", "CONFIRMED", "http://www.boltztrap.org/"),
        ("BoltzTraP2", "CONFIRMED", "https://gitlab.com/sousaw/BoltzTraP2"),
        ("AMSET", "VERIFIED", "UNKNOWN"),
        ("LobsterPy", "LOW_CONF", "UNKNOWN"),
        ("COHP", "LOW_CONF", "UNKNOWN"),
        ("Bader", "VERIFIED", "UNKNOWN"),
        ("DDEC", "UNCERTAIN", "UNKNOWN"),
        ("Critic2", "VERIFIED", "https://github.com/aoterodelaroza/critic2"),
        ("Hirshfeld", "UNCERTAIN", "UNKNOWN"),
        ("FEFF", "VERIFIED", "UNKNOWN"),
        ("xspectra", "LOW_CONF", "UNKNOWN"),
        ("exciting-XS", "LOW_CONF", "UNKNOWN"),
        ("FDMNES", "LOW_CONF", "UNKNOWN"),
        ("CRYSOL", "LOW_CONF", "UNKNOWN"),
        ("XSpectraTools", "UNCERTAIN", "UNKNOWN"),
        ("ezSpectra", "LOW_CONF", "UNKNOWN"),
        ("Libwfa", "LOW_CONF", "UNKNOWN"),
        ("DP", "LOW_CONF", "UNKNOWN"),
        ("Magnon-codes", "LOW_CONF", "UNKNOWN"),
        ("Spirit", "LOW_CONF", "UNKNOWN"),
        ("VAMPIRE", "LOW_CONF", "UNKNOWN"),
        ("TB2J", "LOW_CONF", "UNKNOWN"),
        ("Mumax3", "LOW_CONF", "UNKNOWN"),
        ("McPhase", "LOW_CONF", "UNKNOWN"),
        ("VESTA", "CONFIRMED", "https://jp-minerals.org/vesta/en/"),
        ("XCrySDen", "CONFIRMED", "http://www.xcrysden.org/"),
        ("VMD", "LOW_CONF", "UNKNOWN"),
        ("Avogadro", "LOW_CONF", "UNKNOWN"),
        ("STMng", "LOW_CONF", "UNKNOWN"),
        ("JMol", "LOW_CONF", "UNKNOWN"),
        ("PyMOL", "LOW_CONF", "UNKNOWN"),
        ("OVITO", "LOW_CONF", "UNKNOWN"),
        ("AutoBZ.jl", "LOW_CONF", "UNKNOWN"),
        ("yambopy", "LOW_CONF", "UNKNOWN"),
        ("dbaAutomator", "LOW_CONF", "UNKNOWN"),
        ("gpaw-tools", "LOW_CONF", "UNKNOWN"),
        ("ASE-GUI", "LOW_CONF", "UNKNOWN"),
        ("Nanodcal", "LOW_CONF", "UNKNOWN"),
        ("Transiesta", "LOW_CONF", "UNKNOWN"),
        ("Smeagol", "LOW_CONF", "UNKNOWN"),
        ("MIKA", "LOW_CONF", "UNKNOWN"),
        ("KITE", "LOW_CONF", "UNKNOWN"),
    ],
    
    "Frameworks": [
        ("spglib", "LOW_CONF", "UNKNOWN"),
        ("MatPy", "LOW_CONF", "UNKNOWN"),
        ("AiiDA", "CONFIRMED", "https://www.aiida.net/"),
        ("FireWorks", "CONFIRMED", "https://materialsproject.github.io/fireworks/"),
        ("atomate", "CONFIRMED", "UNKNOWN"),
        ("atomate2", "CONFIRMED", "https://materialsproject.github.io/atomate2/"),
        ("custodian", "CONFIRMED", "https://materialsproject.github.io/custodian/"),
        ("jobflow", "LOW_CONF", "UNKNOWN"),
        ("jobflow-remote", "LOW_CONF", "UNKNOWN"),
        ("Luigi", "LOW_CONF", "UNKNOWN"),
        ("Parsl", "LOW_CONF", "UNKNOWN"),
        ("MyQueue", "LOW_CONF", "UNKNOWN"),
        ("Dask", "LOW_CONF", "UNKNOWN"),
        ("Pyiron", "LOW_CONF", "UNKNOWN"),
        ("AiiDA-VASP", "LOW_CONF", "UNKNOWN"),
        ("AiiDA-QuantumESPRESSO", "LOW_CONF", "UNKNOWN"),
        ("AiiDA-wannier90", "LOW_CONF", "UNKNOWN"),
        ("AiiDA-yambo", "LOW_CONF", "UNKNOWN"),
        ("aiida-fleur", "LOW_CONF", "UNKNOWN"),
        ("AiiDA-plugin-registry", "LOW_CONF", "UNKNOWN"),
        ("Materials-Project", "CONFIRMED", "https://materialsproject.org/"),
        ("AFLOW", "VERIFIED", "UNKNOWN"),
        ("OQMD", "VERIFIED", "UNKNOWN"),
        ("NOMAD", "VERIFIED", "https://nomad-lab.eu/"),
        ("Materials-Cloud", "VERIFIED", "https://materialscloud.org/"),
        ("JARVIS", "LOW_CONF", "UNKNOWN"),
        ("C2DB", "LOW_CONF", "UNKNOWN"),
        ("2DMatPedia", "LOW_CONF", "UNKNOWN"),
        ("pymatgen-db", "LOW_CONF", "UNKNOWN"),
        ("qmpy", "LOW_CONF", "UNKNOWN"),
        ("NCD", "UNCERTAIN", "UNKNOWN"),
        ("ASR", "LOW_CONF", "UNKNOWN"),
        ("pymatgen-analysis", "LOW_CONF", "UNKNOWN"),
        ("matminer", "LOW_CONF", "UNKNOWN"),
        ("MAST", "LOW_CONF", "UNKNOWN"),
        ("Jarvis-Tools", "LOW_CONF", "UNKNOWN"),
        ("Signac", "LOW_CONF", "UNKNOWN"),
    ],
    
    "Niche": [
        ("MPWorks", "LOW_CONF", "UNKNOWN"),
        ("emmet", "LOW_CONF", "UNKNOWN"),
        ("maggma", "LOW_CONF", "UNKNOWN"),
        ("Matbench", "LOW_CONF", "UNKNOWN"),
        ("CatApp", "LOW_CONF", "UNKNOWN"),
        ("CatMAP", "LOW_CONF", "UNKNOWN"),
        ("GASpy", "LOW_CONF", "UNKNOWN"),
        ("AFLOW-ML", "LOW_CONF", "UNKNOWN"),
        ("AFLOW-SYM", "LOW_CONF", "UNKNOWN"),
        ("PyLada", "LOW_CONF", "https://github.com/pylada/pylada"),
        ("Stoner", "LOW_CONF", "UNKNOWN"),
        ("cmpy", "LOW_CONF", "UNKNOWN"),
        ("OSF", "LOW_CONF", "UNKNOWN"),
        ("Zenodo", "LOW_CONF", "UNKNOWN"),
        ("DataVerse", "UNCERTAIN", "UNKNOWN"),
        ("MLIP", "LOW_CONF", "UNKNOWN"),
        ("n2p2", "LOW_CONF", "UNKNOWN"),
        ("SIMPLE-NN", "LOW_CONF", "UNKNOWN"),
        ("AMP", "LOW_CONF", "UNKNOWN"),
        ("SchNetPack", "LOW_CONF", "UNKNOWN"),
        ("MACE", "LOW_CONF", "UNKNOWN"),
        ("NequIP", "LOW_CONF", "UNKNOWN"),
        ("Allegro", "LOW_CONF", "UNKNOWN"),
        ("m3gnet", "LOW_CONF", "UNKNOWN"),
        ("QMCPACK-addons", "LOW_CONF", "UNKNOWN"),
    ],
}

def create_skeleton_md(category, tool_name, confidence, resource):
    """Create skeleton markdown file for a tool"""
    
    template = f"""# {tool_name}

## Official Resources
- Homepage: {resource if resource != "UNKNOWN" else "UNKNOWN - Requires verification"}
- Documentation: UNKNOWN - Requires verification
- Source Repository: UNKNOWN - Requires verification
- License: UNKNOWN - Requires verification

## Overview
**Confidence Level**: {confidence}
**Status**: Documentation pending

[TO BE COMPLETED]

## Theoretical Methods
[TO BE COMPLETED - Requires verification from official sources]

## Capabilities (CRITICAL)
[TO BE COMPLETED - Only verified capabilities from official documentation]

**Sources**: Pending verification

## Inputs & Outputs
**Input formats**: [TO BE COMPLETED]

**Output data types**: [TO BE COMPLETED]

## Interfaces & Ecosystem
[TO BE COMPLETED - Requires verification]

## Limitations & Known Constraints
[TO BE COMPLETED - Requires official documentation review]

## Verification & Sources
**Primary sources**: [TO BE VERIFIED]

**Secondary sources**: [TO BE VERIFIED]

**Confidence**: {confidence}

**Verification status**: ⏸️ PENDING
- Official homepage: {'UNKNOWN' if resource == 'UNKNOWN' else 'TO BE VERIFIED'}
- Documentation: TO BE VERIFIED
- Capabilities: TO BE VERIFIED
"""
    
    return template

def main():
    base_dir = Path("/home/niel/git/Indranil2020.github.io/scientific_tools_consolidated")
    
    created_count = 0
    skipped_count = 0
    
    for category, tools in remaining_tools.items():
        category_dir = base_dir / category
        category_dir.mkdir(exist_ok=True)
        
        for tool_name, confidence, resource in tools:
            # Clean filename
            filename = tool_name.replace("/", "-").replace(" ", "-") + ".md"
            filepath = category_dir / filename
            
            # Skip if already exists
            if filepath.exists():
                skipped_count += 1
                print(f"SKIP: {filepath} (already exists)")
                continue
            
            # Create skeleton
            content = create_skeleton_md(category, tool_name, confidence, resource)
            filepath.write_text(content)
            created_count += 1
            print(f"CREATE: {filepath}")
    
    print(f"\n=== Summary ===")
    print(f"Created: {created_count} files")
    print(f"Skipped: {skipped_count} files (already exist)")
    print(f"Total processed: {created_count + skipped_count}")

if __name__ == "__main__":
    main()
