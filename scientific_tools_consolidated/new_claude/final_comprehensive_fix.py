#!/usr/bin/env python3
"""
COMPREHENSIVE CATALOG FIX v2 - Improved table detection
Fixed logic to properly detect and parse BOTH table formats
"""

import re
from collections import OrderedDict, defaultdict
from pathlib import Path
import glob

def is_valid_code_name(name):
    """Check if string is a valid code name (not citation/number)."""
    name = name.strip()
    if not name or len(name) < 2:
        return False
    if name.isdigit():
        return False
    if name[0].isdigit() and not any(c.isalpha() for c in name):
        return False

    citation_patterns = [
        r'^\d+$',
        r'^\d+-\d+$',
        r'^[A-Z]\d+$',
        r'^\d+\.\d+',
        r'^e\d+$',
    ]

    for pattern in citation_patterns:
        if re.match(pattern, name):
            return False
    return True

def extract_all_codes_both_formats(filepath):
    """Extract ALL codes, handling both table formats."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    structure = OrderedDict()
    seen_codes = {}

    current_cat = None
    current_cat_title = None
    current_subcat = None
    current_subcat_title = None
    current_subsec = None

    lines = content.split('\n')
    i = 0

    while i < len(lines):
        line = lines[i]

        # Main category
        cat_match = re.match(r'^#\s+(\d+)\.\s+(.+)$', line)
        if cat_match:
            current_cat = cat_match.group(1)
            current_cat_title = cat_match.group(2).strip()
            current_subcat = None
            current_subsec = None
            i += 1
            continue

        # Subcategory
        subcat_match = re.match(r'^##\s+(\d+\.\d+)\s+(.+)$', line)
        if subcat_match:
            current_subcat = subcat_match.group(1)
            current_subcat_title = subcat_match.group(2).strip()
            current_subsec = None

            key = (current_cat, current_subcat)
            if key not in structure:
                structure[key] = {
                    'cat_title': current_cat_title,
                    'subcat_title': current_subcat_title,
                    'subsections': OrderedDict(),
                    'codes': []
                }
            i += 1
            continue

        # Subsection
        subsec_match = re.match(r'^###\s+(\d+\.\d+\.\d+)\s+(.+)$', line)
        if subsec_match:
            current_subsec = subsec_match.group(1)
            current_subsec_title = subsec_match.group(2).strip()

            if current_cat and current_subcat:
                key = (current_cat, current_subcat)
                if key in structure:
                    structure[key]['subsections'][current_subsec] = current_subsec_title
            i += 1
            continue

        # Detect table: look for header line (starts with |, contains word characters)
        if current_cat and current_subcat and line.strip().startswith('|') and not line.strip().startswith('|---'):
            # Check if this looks like a table header (contains column names)
            if 'Code' in line or 'License' in line:
                # Found a table header
                i += 1

                # Skip separator line (|---|---|...)
                if i < len(lines) and lines[i].strip().startswith('|---'):
                    i += 1

                # Process code rows until blank line or new section
                while i < len(lines):
                    row = lines[i].strip()

                    # Stop conditions
                    if not row or row.startswith('#') or (row.startswith('|') and 'Code' in row):
                        break

                    # Must be a table row (starts with |, contains **)
                    if not row.startswith('|') or '**' not in row:
                        i += 1
                        continue

                    # Extract all bold code names
                    bold_matches = re.findall(r'\*\*([^*]+)\*\*', row)

                    for code_text in bold_matches:
                        code_name = code_text.strip()

                        if '(' in code_name:
                            code_name = code_name.split('(')[0].strip()

                        if not is_valid_code_name(code_name):
                            continue

                        name_lower = code_name.lower()

                        # First occurrence wins
                        if name_lower not in seen_codes:
                            seen_codes[name_lower] = {
                                'name': code_name,
                                'cat': current_cat,
                                'subcat': current_subcat,
                                'subsec': current_subsec
                            }

                            key = (current_cat, current_subcat)
                            if key in structure:
                                structure[key]['codes'].append({
                                    'name': code_name,
                                    'name_lower': name_lower,
                                    'subsection': current_subsec
                                })

                    i += 1

                continue

        i += 1

    return structure, seen_codes

def parse_all_entries(pattern):
    """Parse ALL entries files for code details."""
    code_details = {}

    for filepath in sorted(glob.glob(pattern)):
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        pattern_re = r'\|\s*(\d+)\s*\|\s*\*\*([^*]+)\*\*\s*\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|'

        for match in re.finditer(pattern_re, content):
            code_id = match.group(1).strip()
            code_name = match.group(2).strip()
            license_info = match.group(3).strip()
            website = match.group(4).strip()
            basis = match.group(5).strip()
            publication = match.group(6).strip()
            specialization = match.group(7).strip()

            if not is_valid_code_name(code_name):
                continue

            name_lower = code_name.lower()

            has_details = license_info and license_info != '-'
            existing = code_details.get(name_lower)

            if not existing or (has_details and existing.get('license') == '-'):
                code_details[name_lower] = {
                    'id': code_id,
                    'name': code_name,
                    'license': license_info if license_info else '-',
                    'website': website if website else '-',
                    'basis': basis if basis else '-',
                    'publication': publication if publication else '-',
                    'specialization': specialization if specialization else '-'
                }

    return code_details

def generate_final_output(structure, seen_codes, code_details, output_path):
    """Generate the final comprehensive catalog."""
    output = []

    # Header
    output.append("# COMPREHENSIVE MATERIALS CODES CATALOG")
    output.append("## Complete Edition - All Codes with Full Details & Proper Categorization")
    output.append("### Fixed: Both Table Formats Included | January 2026")
    output.append("")
    output.append("---")
    output.append("")

    # Statistics
    grand_total = len(seen_codes)
    codes_with_details = sum(1 for name in seen_codes if name in code_details and code_details[name]['license'] != '-')
    codes_missing_details = grand_total - codes_with_details

    output.append(f"**Grand Total Codes**: {grand_total}")
    output.append(f"**Codes with Complete Details**: {codes_with_details} ({100*codes_with_details/grand_total:.1f}%)")
    output.append(f"**Codes Missing Details**: {codes_missing_details} ({100*codes_missing_details/grand_total:.1f}%)")
    output.append("")
    output.append("---")
    output.append("")

    table_header = "| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |"
    table_sep = "|---|---|---|---|---|---|---|---|"

    cat_totals = defaultdict(int)
    current_cat = None
    missing_detail_codes = []

    # Process structure
    for (cat_num, subcat_num), data in structure.items():
        # New main category
        if cat_num != current_cat:
            if current_cat is not None:
                output.append("")
                output.append(f"### **Category {current_cat} Total: {cat_totals[current_cat]} codes**")
                output.append("")

            current_cat = cat_num
            output.append("")
            output.append("=" * 100)
            output.append(f"# {cat_num}. {data['cat_title']}")
            output.append("=" * 100)
            output.append("")

        # Subcategory
        subcat_code_count = len(data['codes'])
        cat_totals[cat_num] += subcat_code_count

        output.append(f"## {subcat_num} {data['subcat_title']}")
        output.append(f"**Subcategory Total: {subcat_code_count} codes**")
        output.append("─" * 80)
        output.append("")

        # Subsections
        for subsec_num, subsec_title in data['subsections'].items():
            output.append(f"### {subsec_num} {subsec_title}")
            output.append("")

        # Codes table
        if data['codes']:
            output.append(table_header)
            output.append(table_sep)

            for idx, code_info in enumerate(data['codes'], 1):
                name_lower = code_info['name_lower']
                code_name = code_info['name']

                details = code_details.get(name_lower)

                if details and details['license'] != '-':
                    row = f"| {idx} | {details['id']} | **{details['name']}** | {details['license']} | {details['website']} | {details['basis']} | {details['publication']} | {details['specialization']} |"
                else:
                    missing_detail_codes.append(code_name)
                    row = f"| {idx} | - | **{code_name}** | - | - | - | - | - |"

                output.append(row)

            output.append("")
        else:
            output.append("*No codes in this subcategory*")
            output.append("")

    # Final category total
    if current_cat:
        output.append("")
        output.append(f"### **Category {current_cat} Total: {cat_totals[current_cat]} codes**")
        output.append("")

    # Summary
    output.append("")
    output.append("=" * 100)
    output.append("# SUMMARY STATISTICS")
    output.append("=" * 100)
    output.append("")
    output.append(f"**Total Codes in Catalog**: {grand_total}")
    output.append(f"**Codes with Complete Details**: {codes_with_details} ({100*codes_with_details/grand_total:.1f}%)")
    output.append(f"**Codes Missing Details**: {codes_missing_details}")
    output.append("")

    output.append("## Category Breakdown")
    output.append("")
    output.append("| Category | Title | Code Count |")
    output.append("|---|---|---|")

    for cat_num in sorted(cat_totals.keys(), key=lambda x: int(x)):
        cat_title = ""
        for (c, s), d in structure.items():
            if c == cat_num:
                cat_title = d['cat_title'][:50]
                break
        output.append(f"| {cat_num} | {cat_title} | {cat_totals[cat_num]} |")

    output.append(f"| **GRAND TOTAL** | | **{grand_total}** |")
    output.append("")

    # Codes missing details section
    if missing_detail_codes:
        output.append("## Codes Missing Details (Not in entries*.md files)")
        output.append(f"**Total Missing**: {len(set(missing_detail_codes))}")
        output.append("")
        output.append("```")
        for code in sorted(set(missing_detail_codes)):
            output.append(f"  - {code}")
        output.append("```")
        output.append("")
        output.append("### Verification Note")
        output.append(f"These {len(set(missing_detail_codes))} codes appear in the comprehensive_code_list.md but do not have matching entries in any of the entries*.md files. This means their full details (license, website, basis set, publication, specialization) are not available in the current documentation.")

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(output))

    return grand_total, codes_with_details, codes_missing_details

def main():
    base = Path("/home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/new_claude")

    print("=" * 70)
    print("COMPREHENSIVE CATALOG FIX v2 - IMPROVED TABLE DETECTION")
    print("=" * 70)

    print("\n[1/3] Extracting codes from comprehensive_code_list.md...")
    print("      Detecting and parsing both table formats...")
    structure, seen_codes = extract_all_codes_both_formats(base / "comprehensive_code_list.md")

    total_in_structure = sum(len(d['codes']) for d in structure.values())
    print(f"      Found {len(structure)} subcategories")
    print(f"      Found {len(seen_codes)} unique valid codes")
    print(f"      Codes in structure: {total_in_structure}")

    print("\n[2/3] Parsing entries*.md files for complete details...")
    code_details = parse_all_entries(str(base / "entries*.md"))
    print(f"      Found details for {len(code_details)} codes")

    print("\n[3/3] Generating comprehensive catalog...")
    output_file = base / "unified_all_codes_VERIFIED_ENHANCED.md"

    total, with_details, missing = generate_final_output(structure, seen_codes, code_details, output_file)

    print(f"\n      Output: {output_file}")
    print(f"      Total codes: {total}")
    print(f"      With details: {with_details} ({100*with_details/total:.1f}%)")
    print(f"      Missing details: {missing}")

    file_size = output_file.stat().st_size / 1024
    print(f"      File size: {file_size:.1f} KB")

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)

if __name__ == "__main__":
    main()
