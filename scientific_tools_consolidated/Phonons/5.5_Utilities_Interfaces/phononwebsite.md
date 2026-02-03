# phononwebsite

## Official Resources
- Homepage: https://henriquemiranda.github.io/phononwebsite/
- Source Repository: https://github.com/henriquemiranda/phononwebsite
- License: MIT License

## Overview
phononwebsite is a web-based tool for visualizing lattice vibrations (phonons) in materials. It provides an interactive 3D visualization of phonon modes, allowing users to explore atomic displacements for different phonon branches and q-points.

**Scientific domain**: Phonon visualization, lattice dynamics education  
**Target user community**: Researchers and students visualizing phonon modes

## Theoretical Methods
- Phonon mode visualization
- Eigenvector animation
- 3D crystal structure display
- Interactive q-point selection
- Mode frequency display

## Capabilities (CRITICAL)
- Interactive 3D phonon visualization
- Atomic displacement animations
- Multiple phonon branches
- Q-point path navigation
- Phonopy band.yaml support
- Web-based interface
- No installation required

## Key Strengths

### Web-Based:
- No installation needed
- Cross-platform
- Shareable visualizations
- Educational tool

### Interactive:
- 3D manipulation
- Mode selection
- Animation controls
- Intuitive interface

## Inputs & Outputs
- **Input formats**:
  - Phonopy band.yaml
  - JSON phonon data
  
- **Output data types**:
  - Interactive visualizations
  - Animations
  - Screenshots

## Interfaces & Ecosystem
- **Phonopy**: Primary data source
- **Web browser**: Interface
- **JavaScript/Three.js**: Rendering


## Advanced Features
- **3D visualization**: Interactive crystal structure display
- **Mode animation**: Atomic displacement animations
- **Q-point navigation**: Interactive path selection
- **Branch selection**: Individual phonon branch visualization
- **Three.js rendering**: High-quality WebGL graphics
- **Shareable links**: Direct URL to specific visualizations

## Performance Characteristics
- Web-based: Runs in browser
- No installation required
- Fast rendering with WebGL
- Cross-platform compatible

## Computational Cost
- Phonopy calculation: External (dominant cost)
- phononwebsite: No computation (visualization only)
- Overall: Zero computational overhead

## Best Practices
- Use converged Phonopy calculations
- Export band.yaml with eigenvectors
- Use modern browser for best performance
- Adjust animation speed for clarity

## Limitations & Known Constraints
- Limited analysis features
- Visualization only
- Requires compatible input
- Browser-dependent

## Application Areas
- Phonon mode visualization
- Teaching and education
- Publication figures
- Mode analysis
- Presentations

## Verification & Sources
**Primary sources**:
1. Website: https://henriquemiranda.github.io/phononwebsite/
2. GitHub: https://github.com/henriquemiranda/phononwebsite

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Web interface: ACCESSIBLE
