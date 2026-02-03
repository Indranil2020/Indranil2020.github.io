#!/usr/bin/env python3
"""
Comprehensive merge script that:
1. Extracts ALL codes from comprehensive_code_list.md with their exact categories
2. Extracts ALL code details from entries*.md files
3. Merges them with proper one-to-one mapping
4. Adds sequential numbering within each subcategory
5. Adds totals at category and subcategory levels
6. Identifies and reports codes with missing details
"""

import re
from collections import defaultdict, OrderedDict
from pathlib import Path
import glob

def parse_comprehensive_list(filepath):
    """Parse comprehensive_code_list.md to extract category structure and code assignments."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Structure: {(cat_num, subcat_num): {'title': ..., 'codes': [...]}}
    structure = OrderedDict()
    current_cat = None
    current_cat_title = None
    current_subcat = None
    current_subcat_title = None
    current_subsec = None
    current_subsec_title = None

    lines = content.split('\n')

    for line in lines:
        # Main category: # 1. TITLE
        cat_match = re.match(r'^#\s+(\d+)\.\s+(.+)$', line)
        if cat_match:
            current_cat = cat_match.group(1)
            current_cat_title = cat_match.group(2).strip()
            current_subcat = None
            current_subsec = None
            continue

        # Subcategory: ## 1.1 Title
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

        # Extract code from table row
        if '|' in line and '**' in line and current_cat and current_subcat:
            # Pattern: | ID | **Name** | ...
            match = re.search(r'\|\s*(\d+)\s*\|\s*\*\*([^*]+)\*\*', line)
            if match:
                code_id = match.group(1).strip()
                code_name = match.group(2).strip()

                key = (current_cat, current_subcat)
                if key in structure:
                    structure[key]['codes'].append({
                        'id': code_id,
                        'name': code_name,
                        'subsection': current_subsec
                    })

    return structure

def parse_all_entries_files(pattern):
    """Parse all entries*.md files to extract code details."""
    code_details = {}

    for filepath in glob.glob(pattern):
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Pattern for table rows with full details
        pattern_re = r'\|\s*(\d+)\s*\|\s*\*\*([^*]+)\*\*\s*\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|'

        for match in re.finditer(pattern_re, content):
            code_id = match.group(1).strip()
            code_name = match.group(2).strip()
            license_info = match.group(3).strip()
            website = match.group(4).strip()
            basis = match.group(5).strip()
            publication = match.group(6).strip()
            specialization = match.group(7).strip()

            # Use code name (lowercase) as primary key
            name_lower = code_name.lower()

            # Store with full details
            if name_lower not in code_details or (license_info != '-' and code_details.get(name_lower, {}).get('license') == '-'):
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

def merge_and_generate(structure, code_details, output_path):
    """Merge structure with details and generate the final file."""
    output = []

    # Header
    output.append("# COMPREHENSIVE MATERIALS CODES CATALOG")
    output.append("## Complete Edition with Proper Category Organization")
    output.append("### Verified & Merged from All Sources | January 2026")
    output.append("")
    output.append("---")
    output.append("")

    # Count totals
    total_codes = 0
    codes_with_details = 0
    codes_missing_details = []

    for key, data in structure.items():
        total_codes += len(data['codes'])

    output.append(f"**Total Codes in Catalog**: {total_codes}")
    output.append("")
    output.append("---")
    output.append("")

    # Table template
    table_header = "| # | ID | Code Name | License | Official Website | Basis Set | Primary Publication with DOI | Specialization |"
    table_sep = "|---|---|---|---|---|---|---|---|"

    # Track category totals
    cat_totals = defaultdict(int)

    # Process each category
    current_cat = None

    for (cat_num, subcat_num), data in structure.items():
        # New main category
        if cat_num != current_cat:
            if current_cat is not None:
                # Print previous category total
                output.append("")
                output.append(f"**Category {current_cat} Total: {cat_totals[current_cat]} codes**")
                output.append("")

            current_cat = cat_num
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

        # Add subsections if any
        for subsec_num, subsec_title in data['subsections'].items():
            output.append(f"### {subsec_num} {subsec_title}")
            output.append("")

        if data['codes']:
            output.append(table_header)
            output.append(table_sep)

            # Number codes sequentially within subcategory
            for idx, code_info in enumerate(data['codes'], 1):
                code_name = code_info['name']
                code_id = code_info['id']
                name_lower = code_name.lower()

                # Try to find details - try exact match first, then variations
                details = None
                search_names = [
                    name_lower,
                    name_lower.replace('+', ''),
                    name_lower.replace('-', ''),
                    name_lower.replace(' ', ''),
                    name_lower.split('/')[0] if '/' in name_lower else None,
                    name_lower.split()[0] if ' ' in name_lower else None,
                ]

                for search_name in search_names:
                    if search_name and search_name in code_details:
                        details = code_details[search_name]
                        break

                if details and details['license'] != '-':
                    codes_with_details += 1
                    row = f"| {idx} | {details['id']} | **{details['name']}** | {details['license']} | {details['website']} | {details['basis']} | {details['publication']} | {details['specialization']} |"
                else:
                    codes_missing_details.append(code_name)
                    row = f"| {idx} | {code_id} | **{code_name}** | - | - | - | - | - |"

                output.append(row)

            output.append("")
        else:
            output.append("*No codes in this subcategory*")
            output.append("")

    # Final category total
    if current_cat:
        output.append("")
        output.append(f"**Category {current_cat} Total: {cat_totals[current_cat]} codes**")
        output.append("")

    # Summary section
    output.append("=" * 100)
    output.append("# SUMMARY STATISTICS")
    output.append("=" * 100)
    output.append("")
    output.append(f"**Total Codes**: {total_codes}")
    output.append(f"**Codes with Full Details**: {codes_with_details} ({100*codes_with_details/total_codes:.1f}%)")
    output.append(f"**Codes Missing Details**: {len(codes_missing_details)}")
    output.append("")

    output.append("### Category Breakdown")
    output.append("| Category | Code Count |")
    output.append("|---|---|")
    for cat_num in sorted(cat_totals.keys(), key=lambda x: int(x)):
        cat_title = ""
        for (c, s), d in structure.items():
            if c == cat_num:
                cat_title = d['cat_title']
                break
        output.append(f"| {cat_num}. {cat_title[:50]} | {cat_totals[cat_num]} |")
    output.append(f"| **TOTAL** | **{total_codes}** |")
    output.append("")

    if codes_missing_details:
        output.append("### Codes Missing Details (need verification)")
        output.append("```")
        for code in sorted(codes_missing_details)[:50]:
            output.append(f"  - {code}")
        if len(codes_missing_details) > 50:
            output.append(f"  ... and {len(codes_missing_details) - 50} more")
        output.append("```")

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(output))

    return total_codes, codes_with_details, len(codes_missing_details)

def main():
    base = Path("/home/niel/git/Indranil2020.github.io/scientific_tools_consolidated/new_claude")

    print("=" * 70)
    print("COMPREHENSIVE CODE CATALOG MERGE")
    print("=" * 70)

    # Step 1: Parse comprehensive list structure
    print("\n[1/3] Parsing comprehensive_code_list.md...")
    structure = parse_comprehensive_list(base / "comprehensive_code_list.md")

    total_in_structure = sum(len(d['codes']) for d in structure.values())
    print(f"      Found {len(structure)} subcategories with {total_in_structure} total code entries")

    # Step 2: Parse all entries files
    print("\n[2/3] Parsing entries*.md files...")
    code_details = parse_all_entries_files(str(base / "entries*.md"))
    print(f"      Found {len(code_details)} unique codes with details")

    # Step 3: Merge and generate
    print("\n[3/3] Merging and generating output...")
    output_file = base / "unified_all_codes_VERIFIED_ENHANCED.md"

    total, with_details, missing = merge_and_generate(structure, code_details, output_file)

    print(f"\n      Output: {output_file}")
    print(f"      Total codes: {total}")
    print(f"      With details: {with_details} ({100*with_details/total:.1f}%)")
    print(f"      Missing details: {missing}")

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)

if __name__ == "__main__":
    main()
