#!/usr/bin/env python3
"""
COMPLETE CATALOG MERGE - Smart extraction from ALL sources
1. Parse ALL entries*.md files as primary source (274 codes with full details)
2. Parse comprehensive_code_list.md for category mapping
3. Merge with proper categorization
4. Add all codes from both sources
"""

import re
from collections import OrderedDict, defaultdict
from pathlib import Path
import glob

def parse_all_entries_comprehensive(pattern):
    """Parse ALL entries files - extract every code with details."""
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

            name_lower = code_name.lower()

            # Store/update with best available details
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

def extract_comprehensive_structure(filepath):
    """Extract category structure from comprehensive list."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    structure = OrderedDict()
    code_to_cat = {}  # Map code names to their primary categories

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

        # Extract codes from table rows - ANY bold text in table
        if '|' in line and '**' in line and current_cat and current_subcat:
            bold_matches = re.findall(r'\*\*([^*]+)\*\*', line)

            for code_text in bold_matches:
                code_name = code_text.strip()

                if not code_name or len(code_name) < 2:
                    continue

                # Remove descriptions in parentheses
                if '(' in code_name:
                    code_name = code_name.split('(')[0].strip()

                # Skip pure numbers
                if code_name.isdigit() or (code_name[0].isdigit() and not any(c.isalpha() for c in code_name)):
                    continue

                # Skip citation patterns
                if re.match(r'^\d+-\d+$|^[A-Z]\d+$|^\d+\.\d+|^e\d+$', code_name):
                    continue

                name_lower = code_name.lower()

                # Store category mapping for this code (first occurrence wins)
                if name_lower not in code_to_cat:
                    code_to_cat[name_lower] = {
                        'cat': current_cat,
                        'subcat': current_subcat,
                        'name': code_name
                    }

                    # Add to structure
                    key = (current_cat, current_subcat)
                    if key in structure:
                        structure[key]['codes'].append({
                            'name': code_name,
                            'name_lower': name_lower,
                            'subsection': current_subsec
                        })

    return structure, code_to_cat

def generate_complete_catalog(structure, code_details, code_to_cat, output_path):
    """Generate complete catalog with ALL codes."""
    output = []

    # Header
    output.append("# COMPREHENSIVE MATERIALS CODES CATALOG")
    output.append("## Complete Edition - All Codes with Full Details")
    output.append("### Comprehensive Merge | January 2026")
    output.append("")
    output.append("---")
    output.append("")

    # Calculate totals
    grand_total = len(code_details)
    codes_with_license = sum(1 for d in code_details.values() if d['license'] != '-')

    output.append(f"**Grand Total Codes**: {grand_total}")
    output.append(f"**Codes with License Info**: {codes_with_license} ({100*codes_with_license/grand_total:.1f}%)")
    output.append("")
    output.append("---")
    output.append("")

    table_header = "| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |"
    table_sep = "|---|---|---|---|---|---|---|---|"

    cat_totals = defaultdict(int)
    current_cat = None

    # Process structure first (categorized codes)
    processed_codes = set()

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
        subcat_total = 0
        table_codes = []

        # Collect codes for this subcategory
        for code_info in data['codes']:
            name_lower = code_info['name_lower']
            code_name = code_info['name']

            details = code_details.get(name_lower)
            if details:
                table_codes.append((code_name, name_lower, details))
                processed_codes.add(name_lower)
                subcat_total += 1

        cat_totals[cat_num] += subcat_total

        output.append(f"## {subcat_num} {data['subcat_title']}")
        output.append(f"**Subcategory Total: {subcat_total} codes**")
        output.append("─" * 80)
        output.append("")

        # Subsections
        for subsec_num, subsec_title in data['subsections'].items():
            output.append(f"### {subsec_num} {subsec_title}")
            output.append("")

        # Code table
        if table_codes:
            output.append(table_header)
            output.append(table_sep)

            for idx, (code_name, name_lower, details) in enumerate(table_codes, 1):
                row = f"| {idx} | {details['id']} | **{details['name']}** | {details['license']} | {details['website']} | {details['basis']} | {details['publication']} | {details['specialization']} |"
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

    # Codes not in comprehensive list (add to category 11)
    uncat_codes = [code for code in code_details.keys() if code not in processed_codes]

    if uncat_codes:
        output.append("")
        output.append("=" * 100)
        output.append("# 12. ADDITIONAL CODES NOT IN COMPREHENSIVE STRUCTURE")
        output.append("=" * 100)
        output.append("")
        output.append("## 12.1 Uncategorized Codes from Entries Files")
        output.append(f"**Subcategory Total: {len(uncat_codes)} codes**")
        output.append("─" * 80)
        output.append("")

        output.append(table_header)
        output.append(table_sep)

        for idx, name_lower in enumerate(sorted(uncat_codes), 1):
            details = code_details[name_lower]
            row = f"| {idx} | {details['id']} | **{details['name']}** | {details['license']} | {details['website']} | {details['basis']} | {details['publication']} | {details['specialization']} |"
            output.append(row)

        output.append("")
        output.append(f"### **Category 12 Total: {len(uncat_codes)} codes**")
        output.append("")

    # Summary
    output.append("")
    output.append("=" * 100)
    output.append("# SUMMARY STATISTICS")
    output.append("=" * 100)
    output.append("")
    output.append(f"**Total Codes**: {grand_total}")
    output.append(f"**Codes with License Information**: {codes_with_license} ({100*codes_with_license/grand_total:.1f}%)")
    output.append("")

    output.append("## Category Breakdown")
    output.append("")
    output.append("| Category | Code Count |")
    output.append("|---|---|")

    for cat_num in sorted(cat_totals.keys(), key=lambda x: int(x)):
        cat_title = ""
        for (c, s), d in structure.items():
            if c == cat_num:
                cat_title = d['cat_title'][:50]
                break
        output.append(f"| {cat_num}. {cat_title} | {cat_totals[cat_num]} |")

    if uncat_codes:
        output.append(f"| 12. Additional Uncategorized Codes | {len(uncat_codes)} |")

    output.append(f"| **GRAND TOTAL** | **{grand_total}** |")
    output.append("")

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(output))

    return grand_total, codes_with_license, len(uncat_codes)

def main():
    base = Path("/home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/new_claude")

    print("=" * 70)
    print("COMPLETE CATALOG MERGE - ALL SOURCES")
    print("=" * 70)

    print("\n[1/3] Parsing ALL entries*.md files...")
    code_details = parse_all_entries_comprehensive(str(base / "entries*.md"))
    print(f"      Found {len(code_details)} unique codes with full details")

    print("\n[2/3] Extracting category structure from comprehensive_code_list.md...")
    structure, code_to_cat = extract_comprehensive_structure(base / "comprehensive_code_list.md")
    print(f"      Found {len(structure)} subcategories")
    print(f"      Mapped {len(code_to_cat)} codes to categories")

    print("\n[3/3] Generating complete merged catalog...")
    output_file = base / "unified_all_codes_VERIFIED_ENHANCED.md"

    total, with_details, uncategorized = generate_complete_catalog(structure, code_details, code_to_cat, output_file)

    print(f"\n      Output: {output_file}")
    print(f"      Total codes: {total}")
    print(f"      With details: {with_details} ({100*with_details/total:.1f}%)")
    print(f"      Uncategorized: {uncategorized}")

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)

if __name__ == "__main__":
    main()
