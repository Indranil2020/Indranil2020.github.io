#!/usr/bin/env python3
"""
FINAL ACCURATE CATALOG - Using comprehensive_code_list.md as source of truth
Extracts ALL 434 codes from comprehensive list with proper categories
Enriches with details from entries*.md files (275 codes have details)
Shows all codes even if details are incomplete
"""

import re
from collections import OrderedDict, defaultdict
from pathlib import Path
import glob

def parse_entries_for_details(pattern):
    """Parse entries files and extract code details."""
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

            name_lower = code_name.lower()

            # Store with best details
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

def extract_comprehensive_structure_with_codes(filepath):
    """Extract COMPLETE structure from comprehensive list including all codes."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    structure = OrderedDict()
    seen_codes = set()  # Track all codes we've seen

    current_cat = None
    current_cat_title = None
    current_subcat = None
    current_subcat_title = None
    current_subsec = None

    lines = content.split('\n')

    for line in lines:
        # Main category
        cat_match = re.match(r'^#\s+(\d+)\.\s+(.+)$', line)
        if cat_match:
            current_cat = cat_match.group(1)
            current_cat_title = cat_match.group(2).strip()
            current_subcat = None
            continue

        # Subcategory
        subcat_match = re.match(r'^##\s+((\d+\.\d+))\s+(.+)$', line)
        if subcat_match:
            current_subcat = subcat_match.group(2)
            current_subcat_title = subcat_match.group(3).strip()
            current_subsec = None

            key = (current_cat, current_subcat)
            if key not in structure:
                structure[key] = {
                    'cat_title': current_cat_title,
                    'subcat_title': current_subcat_title,
                    'subsections': OrderedDict(),
                    'codes': []
                }
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
            continue

        # Extract ALL codes from table rows
        if '|' in line and '**' in line and current_cat and current_subcat:
            bold_matches = re.findall(r'\*\*([^*]+)\*\*', line)

            for code_text in bold_matches:
                code_name = code_text.strip()

                if not code_name or len(code_name) < 2:
                    continue

                if '(' in code_name:
                    code_name = code_name.split('(')[0].strip()

                # Skip pure numbers
                if code_name.isdigit() or (code_name[0].isdigit() and not any(c.isalpha() for c in code_name)):
                    continue

                # Skip citation patterns
                if re.match(r'^\d+-\d+$|^[A-Z]\d+$|^\d+\.\d+|^e\d+$', code_name):
                    continue

                name_lower = code_name.lower()

                # Add EVERY code only once (first occurrence)
                if name_lower not in seen_codes:
                    seen_codes.add(name_lower)

                    key = (current_cat, current_subcat)
                    if key in structure:
                        structure[key]['codes'].append({
                            'name': code_name,
                            'name_lower': name_lower,
                            'subsection': current_subsec
                        })

    return structure, seen_codes

def generate_comprehensive_catalog(structure, seen_codes, code_details, output_path):
    """Generate comprehensive catalog with ALL codes from comprehensive list."""
    output = []

    # Header
    output.append("# COMPREHENSIVE MATERIALS CODES CATALOG")
    output.append("## Complete Edition - All Codes Properly Categorized")
    output.append("### Source: comprehensive_code_list.md | Enriched with entries*.md details | January 2026")
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
    subcat_totals = defaultdict(int)
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
        subcat_key = (cat_num, subcat_num)
        subcat_totals[subcat_key] = subcat_code_count

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
                    # Show with empty details
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
                cat_title = d['cat_title'][:45]
                break
        output.append(f"| {cat_num} | {cat_title} | {cat_totals[cat_num]} |")

    output.append(f"| **GRAND TOTAL** | | **{grand_total}** |")
    output.append("")

    # Missing details list
    if missing_detail_codes:
        output.append("## Codes Missing Details (Not in entries*.md files)")
        output.append(f"**Total Missing**: {len(set(missing_detail_codes))}")
        output.append("")
        output.append("```")
        for code in sorted(set(missing_detail_codes)):
            output.append(f"  - {code}")
        output.append("```")

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(output))

    return grand_total, codes_with_details, codes_missing_details

def main():
    base = Path("/home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/new_claude")

    print("=" * 70)
    print("FINAL ACCURATE CATALOG - USING COMPREHENSIVE LIST AS SOURCE OF TRUTH")
    print("=" * 70)

    print("\n[1/2] Extracting COMPLETE structure from comprehensive_code_list.md...")
    structure, seen_codes = extract_comprehensive_structure_with_codes(base / "comprehensive_code_list.md")
    print(f"      Found {len(structure)} subcategories")
    print(f"      Found {len(seen_codes)} UNIQUE CODES")

    print("\n[2/2] Enriching with details from entries*.md files...")
    code_details = parse_entries_for_details(str(base / "entries*.md"))
    print(f"      Found details for {len(code_details)} codes")

    print("\n[3/3] Generating complete catalog...")
    output_file = base / "unified_all_codes_VERIFIED_ENHANCED.md"
    total, with_details, missing = generate_comprehensive_catalog(structure, seen_codes, code_details, output_file)

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
