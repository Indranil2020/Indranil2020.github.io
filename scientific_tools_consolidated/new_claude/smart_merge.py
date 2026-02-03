#!/usr/bin/env python3
"""
Smart merge script that:
1. Extracts ALL codes from comprehensive_code_list.md with their actual categories
2. Properly handles the document structure (# Main, ## Sub, ### SubSec)
3. Gets details from entries*.md files
4. Merges with proper one-to-one mapping
5. Adds sequential numbering and totals
"""

import re
from collections import OrderedDict, defaultdict
from pathlib import Path
import glob

def extract_all_codes_properly(filepath):
    """Extract all codes with their proper category assignments."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Structure: {(cat_num, subcat_num): {...}}
    structure = OrderedDict()

    current_cat = None
    current_cat_title = None
    current_subcat = None
    current_subcat_title = None
    current_subsec = None
    current_subsec_title = None

    # Track unique codes by lowercase name to avoid duplicates
    seen_codes = {}

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

        # Extract codes from table rows
        # Handle both formats:
        # Format 1: | # | **Code Name** | ... (with numbering)
        # Format 2: | Code | **Code Name** | ... (without numbering)
        if '|' in line and '**' in line and current_cat and current_subcat:
            # Skip header rows and separator rows
            if line.strip().startswith('|---') or 'Code' in line and 'License' in line:
                continue

            # Find ALL bold code names in the line
            bold_matches = re.findall(r'\*\*([^*]+)\*\*', line)

            for code_text in bold_matches:
                code_name = code_text.strip()

                # Skip empty strings and descriptions
                if not code_name or len(code_name) < 2:
                    continue

                # Remove description in parentheses if present
                if '(' in code_name:
                    code_name = code_name.split('(')[0].strip()

                # Skip pure numbers (citations) and patterns
                if code_name.isdigit():
                    continue
                if code_name[0].isdigit() and not any(c.isalpha() for c in code_name):
                    continue

                # Skip citation patterns
                if re.match(r'^\d+-\d+$|^[A-Z]\d+$|^\d+\.\d+|^e\d+$', code_name):
                    continue

                name_lower = code_name.lower()

                # Only add if not seen before (first occurrence in category wins)
                if name_lower not in seen_codes:
                    seen_codes[name_lower] = True

                    key = (current_cat, current_subcat)
                    if key in structure:
                        structure[key]['codes'].append({
                            'name': code_name,
                            'name_lower': name_lower,
                            'subsection': current_subsec
                        })

    return structure, seen_codes

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

def generate_output(structure, code_details, output_path):
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
    grand_total = sum(len(data['codes']) for data in structure.values())
    codes_with_details = 0

    # Count codes with details first pass
    for data in structure.values():
        for code_info in data['codes']:
            name_lower = code_info['name_lower']
            if name_lower in code_details and code_details[name_lower]['license'] != '-':
                codes_with_details += 1

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
        subcat_total = len(data['codes'])
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
                    missing_details.append(code_name)
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
    print("SMART COMPREHENSIVE CODE CATALOG MERGE")
    print("=" * 70)

    print("\n[1/3] Extracting all codes from comprehensive_code_list.md...")
    structure, seen_codes = extract_all_codes_properly(base / "comprehensive_code_list.md")

    total_codes_in_structure = sum(len(d['codes']) for d in structure.values())
    print(f"      Found {len(structure)} subcategories")
    print(f"      Found {total_codes_in_structure} unique valid codes")

    print("\n[2/3] Parsing entries*.md files...")
    code_details = parse_all_entries(str(base / "entries*.md"))
    print(f"      Found details for {len(code_details)} codes")

    print("\n[3/3] Generating output...")
    output_file = base / "unified_all_codes_VERIFIED_ENHANCED.md"

    total, with_details = generate_output(structure, code_details, output_file)

    print(f"\n      Output: {output_file}")
    print(f"      Total codes: {total}")
    print(f"      With details: {with_details} ({100*with_details/total:.1f}%)")
    print(f"      Missing details: {total - with_details}")

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)

if __name__ == "__main__":
    main()
