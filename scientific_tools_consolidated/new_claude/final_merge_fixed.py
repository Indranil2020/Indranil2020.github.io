#!/usr/bin/env python3
"""
Fixed merge script that properly extracts code names (not citation numbers).
"""

import re
from collections import OrderedDict, defaultdict
from pathlib import Path
import glob

def is_valid_code_name(name):
    """Check if a string is likely a valid code name (not a citation number)."""
    name = name.strip()

    # Reject if empty or too short
    if not name or len(name) < 2:
        return False

    # Reject if it's just a number
    if name.isdigit():
        return False

    # Reject if it starts with a number and has no letters
    if name[0].isdigit() and not any(c.isalpha() for c in name):
        return False

    # Reject common citation patterns
    citation_patterns = [
        r'^\d+$',           # Pure number
        r'^\d+-\d+$',       # Range like 558-561
        r'^[A-Z]\d+$',      # Like B47
        r'^\d+\.\d+',       # Like 47.558
        r'^e\d+$',          # Like e1606
    ]

    for pattern in citation_patterns:
        if re.match(pattern, name):
            return False

    return True

def extract_codes_from_comprehensive(filepath):
    """Extract all unique codes with their primary categories."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    current_cat = None
    current_cat_title = None
    current_subcat = None
    current_subcat_title = None
    current_subsec = None
    current_subsec_title = None

    code_categories = {}
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

        # Extract codes from table rows - match the actual format: | # | **Code** (description) |
        if '|' in line and current_cat and current_subcat:
            # Pattern: | # | **Code Name** possibly with (description) | ...
            # The code name is in bold in the second column
            match = re.match(r'\|\s*(\d+)\s*\|\s*\*\*([^*\(]+)\*\*', line)
            if match:
                code_id = match.group(1).strip()
                code_name = match.group(2).strip()

                # Clean up code name (remove descriptions in parentheses)
                if '(' in code_name:
                    code_name = code_name.split('(')[0].strip()

                if not is_valid_code_name(code_name):
                    continue

                name_lower = code_name.lower()

                # Only assign to first category found
                if name_lower not in code_categories:
                    code_categories[name_lower] = {
                        'name': code_name,
                        'id': code_id,
                        'category': current_cat,
                        'subcategory': current_subcat
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

        # Pattern for table rows
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

def generate_output(structure, code_categories, code_details, output_path):
    """Generate the final merged file."""
    output = []

    # Header
    output.append("# COMPREHENSIVE MATERIALS CODES CATALOG")
    output.append("## Complete Edition - All Codes Organized with Numbering & Totals")
    output.append("### Verified & Merged | January 2026")
    output.append("")
    output.append("---")
    output.append("")

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

    cat_totals = defaultdict(int)
    missing_details = []
    current_cat = None

    for (cat_num, subcat_num), data in structure.items():
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

        subcat_total = len(data['codes'])
        cat_totals[cat_num] += subcat_total

        output.append(f"## {subcat_num} {data['subcat_title']}")
        output.append(f"**Subcategory Total: {subcat_total} codes**")
        output.append("─" * 80)
        output.append("")

        for subsec_num, subsec_title in data['subsections'].items():
            output.append(f"### {subsec_num} {subsec_title}")
            output.append("")

        if data['codes']:
            output.append(table_header)
            output.append(table_sep)

            for idx, code_info in enumerate(data['codes'], 1):
                name_lower = code_info['name_lower']
                code_name = code_info['name']
                orig_id = code_info['id']

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
    print("FIXED COMPREHENSIVE CODE CATALOG MERGE")
    print("=" * 70)

    print("\n[1/3] Extracting codes from comprehensive_code_list.md...")
    structure, code_categories = extract_codes_from_comprehensive(base / "comprehensive_code_list.md")

    print(f"      Found {len(structure)} subcategories")
    print(f"      Found {len(code_categories)} unique valid codes")

    print("\n[2/3] Parsing entries*.md files...")
    code_details = parse_all_entries(str(base / "entries*.md"))
    print(f"      Found details for {len(code_details)} codes")

    print("\n[3/3] Generating output...")
    output_file = base / "unified_all_codes_VERIFIED_ENHANCED.md"

    total, with_details = generate_output(structure, code_categories, code_details, output_file)

    print(f"\n      Output: {output_file}")
    print(f"      Total codes: {total}")
    print(f"      With details: {with_details} ({100*with_details/total:.1f}%)")
    print(f"      Missing: {total - with_details}")

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)

if __name__ == "__main__":
    main()
