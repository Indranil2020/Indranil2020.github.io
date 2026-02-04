# findmagsym

## Official Resources
- GitHub: https://github.com/yuanlinding/findmagsym
- Web Application: Available online
- License: Check repository

## Overview
findmagsym is a web-based application for finding the magnetic space group of a magnetic crystal structure. It takes a crystal structure with magnetic moments as input and determines the appropriate magnetic space group, generating MCIF (Magnetic CIF) format output files compatible with visualization and analysis tools.

**Scientific domain**: Magnetic crystallography, magnetic space groups, structure analysis
**Target user community**: Researchers studying magnetic materials and needing magnetic symmetry classification

## Theoretical Methods
- Magnetic space group determination
- Magnetic symmetry operation analysis
- Time-reversal symmetry combinations
- Magnetic moment transformation rules
- BNS and OG notation support

## Capabilities (CRITICAL)
- **MSG Determination**: Automatic magnetic space group finding
- **MCIF Output**: Magnetic CIF file generation
- **Web Interface**: Browser-based accessibility
- **Multiple Notations**: BNS and OG conventions
- **Visualization**: Structure visualization support
- **Moment Analysis**: Magnetic moment symmetry checking

**Sources**: GitHub repository, web application

## Key Strengths

### Web-Based:
- No installation required
- Browser accessible
- Immediate results
- User-friendly interface

### MCIF Support:
- Standard output format
- Compatible with VESTA
- Bilbao server compatible
- Documentation standard

### Practical Focus:
- Real structure input
- Moment handling
- Multiple conventions
- Error checking

## Inputs & Outputs
- **Input formats**:
  - CIF files with moments
  - Structure + moment vectors
  - Manual input option
  
- **Output data types**:
  - Magnetic space group number
  - MCIF format file
  - Symmetry operations
  - Magnetic point group

## Installation
Web-based: No installation needed, access via browser.

For local use:
```bash
git clone https://github.com/yuanlinding/findmagsym.git
cd findmagsym
# Follow repository instructions
```

## Usage Examples
1. Navigate to web application
2. Upload CIF file with magnetic structure
3. Specify magnetic moments
4. Click "Find Magnetic Space Group"
5. Download MCIF output

## Performance Characteristics
- **Speed**: Fast symmetry analysis
- **Accuracy**: Validated against known structures
- **Accessibility**: Web-based, no setup

## Limitations & Known Constraints
- **Web dependency**: Requires internet for web version
- **Input format**: Specific format requirements
- **Complex structures**: May need manual verification

## Comparison with Other Tools
- **vs FINDSYM**: findmagsym magnetic-focused
- **vs Bilbao MAGNDATA**: Different interface approach
- **vs STRCONVERT**: findmagsym web-based
- **Unique strength**: Simple web interface for MSG determination

## Application Areas
- Magnetic structure refinement
- Neutron diffraction analysis
- Magnetic phase identification
- Symmetry-based property prediction
- Database entry preparation

## Best Practices
- Verify input structure accuracy
- Check moment directions carefully
- Compare BNS and OG notations
- Validate against known materials

## Community and Support
- GitHub repository
- Web application support
- Academic development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/yuanlinding/findmagsym

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Web application: AVAILABLE
- Method: Magnetic space group determination
