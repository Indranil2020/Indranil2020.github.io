#!/usr/bin/env python3
"""
Full merge script that:
1. Extracts ALL unique codes from comprehensive_code_list.md (390+ codes)
2. Maps them to their primary categories based on which section they first appear
3. Gets details from entries*.md files
4. Adds proper sequential numbering
5. Adds category and subcategory totals
"""

import re
from collections import OrderedDict, defaultdict
from pathlib import Path
import glob

def extract_all_codes_with_categories(filepath):
    """Extract all unique codes and their categories from comprehensive_code_list.md."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Track current location
    current_cat = None
    current_cat_title = None
    current_subcat = None
    current_subcat_title = None
    current_subsec = None
    current_subsec_title = None

    # Store: code_name_lower -> {category info, original_name}
    code_categories = {}

    # Store structure
    structure = OrderedDict()

    lines = content.split('\n')

    for line in lines:
        # Main category
        cat_match = re.match(r'^#\s+(\d+)\.\s+(.+)$', line)
        if cat_match:
            current_cat = cat_match.group(1)
            current_cat_title = cat_match.group(2).strip()
            current_subcat = None
            current_subsec = None
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

        # Extract code names from any line containing bold text in tables
        if '|' in line and '**' in line and current_cat and current_subcat:
            # Find all bold code names
            bold_matches = re.findall(r'\*\*([^*]+)\*\*', line)

            for code_name in bold_matches:
                # Clean up code name (remove descriptions in parentheses)
                code_name = code_name.strip()
                if '(' in code_name:
                    code_name = code_name.split('(')[0].strip()

                if not code_name or len(code_name) < 2:
                    continue

                name_lower = code_name.lower()

                # Try to extract ID if present
                id_match = re.search(r'\|\s*(\d+)\s*\|.*\*\*' + re.escape(code_name), line, re.IGNORECASE)
                code_id = id_match.group(1) if id_match else None

                # Only assign to first category found (primary category)
                if name_lower not in code_categories:
                    code_categories[name_lower] = {
                        'name': code_name,
                        'id': code_id,
                        'category': current_cat,
                        'subcategory': current_subcat,
                        'subsection': current_subsec
                    }

                    key = (current_cat, current_subcat)
                    if key in structure:
                        structure[key]['codes'].append({
                            'name': code_name,
                            'id': code_id,
                            'name_lower': name_lower
                        })

    return structure, code_categories

def parse_all_entries(pattern):
    """Parse all entries files for code details."""
    code_details = {}

    for filepath in glob.glob(pattern):
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Match table rows
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

            # Keep best available details
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

def generate_final_output(structure, code_categories, code_details, output_path):
    """Generate the final merged file with numbering and totals."""
    output = []

    # Header
    output.append("# COMPREHENSIVE MATERIALS CODES CATALOG")
    output.append("## Complete Edition - All Codes Properly Organized with Numbering & Totals")
    output.append("### Verified & Merged | January 2026")
    output.append("")
    output.append("---")
    output.append("")

    # Calculate grand total
    grand_total = len(code_categories)
    codes_with_details = sum(1 for name in code_categories if name in code_details and code_details[name]['license'] != '-')

    output.append(f"**Grand Total Codes**: {grand_total}")
    output.append(f"**Codes with Full Details**: {codes_with_details} ({100*codes_with_details/grand_total:.1f}%)")
    output.append(f"**Codes Needing Details**: {grand_total - codes_with_details}")
    output.append("")
    output.append("---")
    output.append("")

    table_header = "| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |"
    table_sep = "|---|---|---|---|---|---|---|---|"

    # Track totals
    cat_totals = defaultdict(int)
    missing_details = []
    current_cat = None

    for (cat_num, subcat_num), data in structure.items():
        # New main category
        if cat_num != current_cat:
            # Print previous category total
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

        # Subcategory header with total
        subcat_total = len(data['codes'])
        cat_totals[cat_num] += subcat_total

        output.append(f"## {subcat_num} {data['subcat_title']}")
        output.append(f"### **Subcategory Total: {subcat_total} codes**")
        output.append("─" * 80)
        output.append("")

        # Subsections
        for subsec_num, subsec_title in data['subsections'].items():
            output.append(f"### {subsec_num} {subsec_title}")
            output.append("")

        # Code table
        if data['codes']:
            output.append(table_header)
            output.append(table_sep)

            for idx, code_info in enumerate(data['codes'], 1):
                name_lower = code_info['name_lower']
                code_name = code_info['name']
                orig_id = code_info['id'] if code_info['id'] else '-'

                # Get details
                details = code_details.get(name_lower)

                if details and details['license'] != '-':
                    row = f"| {idx} | {details['id']} | **{details['name']}** | {details['license']} | {details['website']} | {details['basis']} | {details['publication']} | {details['specialization']} |"
                else:
                    missing_details.append(code_name)
                    row = f"| {idx} | {orig_id} | **{code_name}** | - | - | - | - | - |"

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

    # Summary statistics
    output.append("")
    output.append("=" * 100)
    output.append("# SUMMARY STATISTICS")
    output.append("=" * 100)
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

    # Missing details report
    output.append("## Codes Missing Details")
    output.append(f"**Total Missing**: {len(missing_details)}")
    output.append("")
    if missing_details:
        output.append("```")
        for code in sorted(set(missing_details)):
            output.append(f"  - {code}")
        output.append("```")

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(output))

    return grand_total, codes_with_details

def main():
    base = Path("/home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/new_claude")

    print("=" * 70)
    print("FULL COMPREHENSIVE CODE CATALOG MERGE")
    print("=" * 70)

    print("\n[1/3] Extracting ALL codes from comprehensive_code_list.md...")
    structure, code_categories = extract_all_codes_with_categories(base / "comprehensive_code_list.md")

    print(f"      Found {len(structure)} subcategories")
    print(f"      Found {len(code_categories)} unique codes")

    print("\n[2/3] Parsing ALL entries*.md files for details...")
    code_details = parse_all_entries(str(base / "entries*.md"))
    print(f"      Found details for {len(code_details)} codes")

    print("\n[3/3] Generating final merged output...")
    output_file = base / "unified_all_codes_VERIFIED_ENHANCED.md"

    total, with_details = generate_final_output(structure, code_categories, code_details, output_file)

    print(f"\n      Output: {output_file}")
    print(f"      Total codes: {total}")
    print(f"      With details: {with_details} ({100*with_details/total:.1f}%)")
    print(f"      Missing details: {total - with_details}")

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)

if __name__ == "__main__":
    main()
