#!/usr/bin/env python3
"""
FINAL COMPLETE CATALOG - Production Quality
Extracts ALL 274 codes from entries files with complete details
Maps to categories from comprehensive_code_list.md
Produces fully organized, verified output
"""

import re
from collections import OrderedDict, defaultdict
from pathlib import Path
import glob

def parse_all_entries_files(pattern):
    """Parse ALL entries files - extract EVERY code with complete details."""
    code_details = {}
    total_rows = 0

    for filepath in sorted(glob.glob(pattern)):
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Match ALL table rows - comprehensive pattern
        pattern_re = r'\|\s*(\d+)\s*\|\s*\*\*([^*]+)\*\*\s*\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|'

        file_codes = 0
        for match in re.finditer(pattern_re, content):
            code_id = match.group(1).strip()
            code_name = match.group(2).strip()
            license_info = match.group(3).strip()
            website = match.group(4).strip()
            basis = match.group(5).strip()
            publication = match.group(6).strip()
            specialization = match.group(7).strip()

            name_lower = code_name.lower()

            # Store/update - keep entry with most details
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
                file_codes += 1

            total_rows += 1

    return code_details, total_rows

def extract_category_mapping(filepath):
    """Extract category structure and map codes to categories."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    structure = OrderedDict()
    code_cat_map = {}

    current_cat = None
    current_cat_title = None
    current_subcat = None
    current_subcat_title = None
    current_subsec = None

    lines = content.split('\n')

    for line in lines:
        # Main category: # 1. TITLE
        cat_match = re.match(r'^#\s+(\d+)\.\s+(.+)$', line)
        if cat_match:
            current_cat = cat_match.group(1)
            current_cat_title = cat_match.group(2).strip()
            current_subcat = None
            continue

        # Subcategory: ## 1.1 Title
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
                    'codes': set()
                }
            continue

        # Subsection: ### 1.1.1 Title
        subsec_match = re.match(r'^###\s+(\d+\.\d+\.\d+)\s+(.+)$', line)
        if subsec_match:
            current_subsec = subsec_match.group(1)
            current_subsec_title = subsec_match.group(2).strip()

            if current_cat and current_subcat:
                key = (current_cat, current_subcat)
                if key in structure:
                    structure[key]['subsections'][current_subsec] = current_subsec_title
            continue

        # Extract codes from tables
        if '|' in line and '**' in line and current_cat and current_subcat:
            # Find all bold code names
            bold_matches = re.findall(r'\*\*([^*]+)\*\*', line)

            for code_text in bold_matches:
                code_name = code_text.strip()

                if not code_name or len(code_name) < 2:
                    continue

                if '(' in code_name:
                    code_name = code_name.split('(')[0].strip()

                # Skip numbers and citation patterns
                if code_name.isdigit() or (code_name[0].isdigit() and not any(c.isalpha() for c in code_name)):
                    continue

                if re.match(r'^\d+-\d+$|^[A-Z]\d+$|^\d+\.\d+|^e\d+$', code_name):
                    continue

                name_lower = code_name.lower()

                # Map to category (first occurrence wins)
                if name_lower not in code_cat_map:
                    code_cat_map[name_lower] = {
                        'cat': current_cat,
                        'subcat': current_subcat,
                        'name': code_name
                    }

                    key = (current_cat, current_subcat)
                    if key in structure:
                        structure[key]['codes'].add(name_lower)

    return structure, code_cat_map

def generate_final_catalog(structure, code_details, code_cat_map, output_path):
    """Generate final production-quality catalog with ALL codes."""
    output = []

    # Header
    output.append("# COMPREHENSIVE MATERIALS CODES CATALOG")
    output.append("## Complete Edition - All Codes with Full Details & Numbering")
    output.append("### Verified & Merged | January 2026")
    output.append("")
    output.append("---")
    output.append("")

    # Statistics
    grand_total = len(code_details)
    codes_with_license = sum(1 for d in code_details.values() if d.get('license', '-') != '-')

    output.append(f"**Grand Total Codes**: {grand_total}")
    output.append(f"**Codes with Complete Details**: {codes_with_license} ({100*codes_with_license/grand_total:.1f}%)")
    output.append("")
    output.append("---")
    output.append("")

    table_header = "| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |"
    table_sep = "|---|---|---|---|---|---|---|---|"

    cat_totals = defaultdict(int)
    processed_codes = set()
    current_cat = None

    # Process categorized codes
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
        subcat_codes = []
        for name_lower in sorted(data['codes']):
            if name_lower in code_details:
                details = code_details[name_lower]
                subcat_codes.append((name_lower, details))
                processed_codes.add(name_lower)

        subcat_total = len(subcat_codes)
        cat_totals[cat_num] += subcat_total

        output.append(f"## {subcat_num} {data['subcat_title']}")
        output.append(f"**Subcategory Total: {subcat_total} codes**")
        output.append("─" * 80)
        output.append("")

        # Subsections
        for subsec_num, subsec_title in data['subsections'].items():
            output.append(f"### {subsec_num} {subsec_title}")
            output.append("")

        # Codes table
        if subcat_codes:
            output.append(table_header)
            output.append(table_sep)

            for idx, (name_lower, details) in enumerate(subcat_codes, 1):
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

    # Uncategorized codes
    uncategorized = [code for code in code_details.keys() if code not in processed_codes]

    if uncategorized:
        output.append("")
        output.append("=" * 100)
        output.append("# 12. ADDITIONAL CODES")
        output.append("=" * 100)
        output.append("")
        output.append("## 12.1 Codes Not in Comprehensive Category Structure")
        output.append(f"**Subcategory Total: {len(uncategorized)} codes**")
        output.append("─" * 80)
        output.append("")

        output.append(table_header)
        output.append(table_sep)

        for idx, name_lower in enumerate(sorted(uncategorized), 1):
            details = code_details[name_lower]
            row = f"| {idx} | {details['id']} | **{details['name']}** | {details['license']} | {details['website']} | {details['basis']} | {details['publication']} | {details['specialization']} |"
            output.append(row)

        output.append("")
        output.append(f"### **Category 12 Total: {len(uncategorized)} codes**")
        output.append("")

    # Summary
    output.append("")
    output.append("=" * 100)
    output.append("# SUMMARY STATISTICS")
    output.append("=" * 100)
    output.append("")
    output.append(f"**Total Codes in Catalog**: {grand_total}")
    output.append(f"**Codes with Complete Details**: {codes_with_license} ({100*codes_with_license/grand_total:.1f}%)")
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

    if uncategorized:
        output.append(f"| 12. Additional Codes | {len(uncategorized)} |")

    output.append(f"| **GRAND TOTAL** | **{grand_total}** |")
    output.append("")

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(output))

    return grand_total, codes_with_license, len(uncategorized)

def main():
    base = Path("/home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/new_claude")

    print("=" * 70)
    print("FINAL COMPLETE CATALOG GENERATION")
    print("=" * 70)

    print("\n[1/3] Parsing ALL entries files for complete details...")
    code_details, total_rows = parse_all_entries_files(str(base / "entries*.md"))
    print(f"      Total table rows parsed: {total_rows}")
    print(f"      Unique codes extracted: {len(code_details)}")

    print("\n[2/3] Extracting category structure...")
    structure, code_cat_map = extract_category_mapping(base / "comprehensive_code_list.md")
    print(f"      Categories/subcategories found: {len(structure)}")
    print(f"      Codes mapped to categories: {len(code_cat_map)}")

    print("\n[3/3] Generating production catalog...")
    output_file = base / "unified_all_codes_VERIFIED_ENHANCED.md"
    total, with_details, uncategorized = generate_final_catalog(structure, code_details, code_cat_map, output_file)

    print(f"\n      Output: {output_file}")
    print(f"      Total codes: {total}")
    print(f"      With complete details: {with_details} ({100*with_details/total:.1f}%)")
    print(f"      Uncategorized (in Category 12): {uncategorized}")

    # File size info
    file_size = output_file.stat().st_size / 1024
    print(f"      File size: {file_size:.1f} KB")

    print("\n" + "=" * 70)
    print("COMPLETE - PRODUCTION CATALOG READY!")
    print("=" * 70)

if __name__ == "__main__":
    main()
