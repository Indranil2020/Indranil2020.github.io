#!/usr/bin/env python3
"""
Merge newly-downloaded PDFs from
  scientific_tools_consolidated/Papers_of_Codes/New_downloaded_Papers/Papers_of_Codes/<Cat>/<Code>/*.pdf
into the live tree at
  scientific_tools_consolidated/Papers_of_Codes/<Cat>/<Code>/*.pdf
and update the corresponding `- Paper:` lines (PLACEHOLDER) inside
  scientific_tools_consolidated/MASTER_LIST_COMPLETE_372.md
so the website reflects the new papers.

Reads:
  - The new-PDFs directory tree (source of truth for filenames).
  - The current master list (matches by code name).
Writes:
  - Existing PDF files moved into Papers_of_Codes/<Cat>/<Code>/
  - Master list entries with PLACEHOLDER lines replaced by real paper links.
  - A migration report MIGRATION_REPORT.md.

Idempotent: safe to re-run; already-merged PDFs are skipped, already-updated
master list lines are not re-touched.
"""
import json
import re
import shutil
from pathlib import Path
from collections import defaultdict
from urllib.parse import unquote, quote

ROOT = Path(__file__).resolve().parent.parent
TOOLS_BASE = ROOT / "scientific_tools_consolidated"
NEW_DIR = TOOLS_BASE / "Papers_of_Codes" / "New_downloaded_Papers" / "Papers_of_Codes"
DEST_BASE = TOOLS_BASE / "Papers_of_Codes"
MASTER = TOOLS_BASE / "MASTER_LIST_COMPLETE_372.md"
REPORT = TOOLS_BASE / "MIGRATION_REPORT.md"

# Master list entry parsing — must match build_tools_db.py's ENTRY_RE exactly so
# every code in the live database is considered.
ENTRY_RE = re.compile(r'^\*\*([0-9]+[a-z]?)\.\s+(.+?)\*\*\*?(?:\s*\(.*?\))?\s*$')
PAPER_LINE_RE = re.compile(r'^- Paper:\s*(.*)$')


def discover_new_pdfs():
    """Return dict: {code_name_normalized: [(cat, code_dir, filename, abs_src), ...]}"""
    out = defaultdict(list)
    for pdf in sorted(NEW_DIR.rglob('*.pdf')):
        rel = pdf.relative_to(NEW_DIR)
        parts = rel.parts
        if len(parts) < 3:
            continue  # need Cat/Code/file.pdf
        cat = parts[0]
        code_dir = parts[1]
        filename = parts[-1]
        out[code_dir.lower()].append({
            'cat': cat,
            'code_dir': code_dir,
            'filename': filename,
            'src': pdf,
        })
    return out


def parse_master_entries(text):
    """Parse master list and return list of dicts:
       {idx, name, line_start, line_end (exclusive), header, paper_line_no (or None), paper_line}.
       line numbers are 0-indexed against the splitlines() array.
    """
    lines = text.splitlines()
    entries = []
    i = 0
    cur = None
    while i < len(lines):
        m = ENTRY_RE.match(lines[i])
        if m:
            if cur is not None:
                cur['line_end'] = i
                entries.append(cur)
            cur = {
                'idx': m.group(1),
                'name': m.group(2).strip(),
                'line_start': i,
                'paper_line_no': None,
                'paper_line': None,
            }
        elif cur is not None:
            pm = PAPER_LINE_RE.match(lines[i])
            if pm and cur['paper_line_no'] is None:
                cur['paper_line_no'] = i
                cur['paper_line'] = lines[i]
            # Also stop the entry when we hit a category header / hr
            if lines[i].startswith('## ') or lines[i].startswith('### ') or lines[i].startswith('---'):
                cur['line_end'] = i
                entries.append(cur)
                cur = None
        i += 1
    if cur is not None:
        cur['line_end'] = len(lines)
        entries.append(cur)
    return entries, lines


def normalize(s):
    s = unquote(str(s or ''))
    return re.sub(r'[^a-z0-9]+', '', s.lower())


def main():
    new_pdfs = discover_new_pdfs()
    print(f"Discovered {sum(len(v) for v in new_pdfs.values())} new PDFs across {len(new_pdfs)} code dirs")

    text = MASTER.read_text()
    entries, lines = parse_master_entries(text)
    print(f"Parsed {len(entries)} entries from master list")

    # Index entries by normalized name for fuzzy matching
    by_norm = defaultdict(list)
    for e in entries:
        by_norm[normalize(e['name'])].append(e)

    matched = []      # (entry, new_pdf_records)
    unmatched_pdfs = []
    moved_count = 0
    already_present_count = 0
    master_updates = 0
    placeholder_replaced = 0

    for code_lower, recs in new_pdfs.items():
        # Try exact normalized match first
        key = normalize(recs[0]['code_dir'])
        candidates = by_norm.get(key, [])
        if not candidates:
            # Try a few common variants: spaces, hyphens, dots stripped already; try without '.jl' suffix
            alt = key.replace('jl', '')
            candidates = by_norm.get(alt, [])
        if not candidates:
            unmatched_pdfs.append((recs[0]['code_dir'], [r['filename'] for r in recs]))
            continue
        # If multiple, prefer the one whose category matches
        chosen = candidates[0]
        if len(candidates) > 1:
            for c in candidates:
                if c.get('paper_line') and recs[0]['cat'].lower() in c['paper_line'].lower():
                    chosen = c
                    break
        matched.append((chosen, recs))

    print(f"Matched {len(matched)} code dirs to master list entries; {len(unmatched_pdfs)} unmatched")

    # Move PDFs and update master list
    for entry, recs in matched:
        cat = recs[0]['cat']
        # URL-decode dir name so e.g. 'Qbox_%28TDDFT' lands as 'Qbox_(TDDFT' on disk;
        # close any unmatched opening paren so dirs are well-formed.
        code_dir = unquote(recs[0]['code_dir'])
        if code_dir.count('(') > code_dir.count(')'):
            code_dir = code_dir + ')' * (code_dir.count('(') - code_dir.count(')'))
        dest_dir = DEST_BASE / cat / code_dir
        dest_dir.mkdir(parents=True, exist_ok=True)

        new_paper_paths = []
        for r in recs:
            dest = dest_dir / unquote(r['filename'])
            if dest.exists() and dest.stat().st_size == r['src'].stat().st_size:
                already_present_count += 1
            else:
                shutil.copy2(r['src'], dest)
                moved_count += 1
            rel_path = dest.relative_to(TOOLS_BASE).as_posix()
            # URL-encode only special chars that break markdown link parsing; keep '/' as is
            rel_path_md = quote(rel_path, safe='/._-')
            new_paper_paths.append((unquote(r['filename']), rel_path_md))

        # Build new Paper line
        link_strs = [f'[{name}]({path})' for name, path in new_paper_paths]
        new_line = '- Paper: ' + ', '.join(link_strs)

        if entry['paper_line_no'] is not None:
            old = lines[entry['paper_line_no']]
            if 'PLACEHOLDER' in old:
                lines[entry['paper_line_no']] = new_line
                master_updates += 1
                placeholder_replaced += 1
            elif old.strip() != new_line.strip():
                # Existing real paper line — only append new ones not already present
                existing_paths = set(re.findall(r'\(([^)]+\.pdf)\)', old))
                additions = [(n, p) for n, p in new_paper_paths if p not in existing_paths]
                if additions:
                    extra = ', '.join(f'[{n}]({p})' for n, p in additions)
                    lines[entry['paper_line_no']] = old.rstrip() + ', ' + extra
                    master_updates += 1
        else:
            # No Paper line yet: insert one right before line_end
            insert_at = entry['line_end']
            lines.insert(insert_at, new_line)
            master_updates += 1
            # Adjust subsequent line numbers in entries
            for ent2 in entries:
                if ent2['line_start'] > insert_at:
                    ent2['line_start'] += 1
                if ent2['paper_line_no'] is not None and ent2['paper_line_no'] >= insert_at:
                    ent2['paper_line_no'] += 1
                ent2['line_end'] = ent2.get('line_end', 0) + (1 if ent2['line_end'] >= insert_at else 0)

    # Write master list back
    new_text = '\n'.join(lines)
    if not new_text.endswith('\n'):
        new_text += '\n'
    if new_text != text:
        MASTER.write_text(new_text)
        print(f"Updated master list ({master_updates} lines, {placeholder_replaced} PLACEHOLDERs replaced)")
    else:
        print("Master list already up to date")

    # Re-parse to compute final stats
    entries2, _ = parse_master_entries(MASTER.read_text())
    total = len(entries2)
    placeholders = sum(1 for e in entries2 if e['paper_line'] and 'PLACEHOLDER' in e['paper_line'])
    no_paper_line = sum(1 for e in entries2 if e['paper_line_no'] is None)
    with_real = total - placeholders - no_paper_line

    # Write report
    rep_lines = [
        "# Paper Migration Report",
        "",
        f"Generated by `scripts/merge_new_papers.py`.",
        "",
        "## Migration outcome",
        "",
        f"- New PDFs discovered:        **{sum(len(v) for v in new_pdfs.values())}**",
        f"- Code dirs matched to master entries: **{len(matched)}**",
        f"- Code dirs unmatched:        **{len(unmatched_pdfs)}**",
        f"- PDFs copied into live tree: **{moved_count}**",
        f"- PDFs already present:       **{already_present_count}**",
        f"- Master list lines updated:  **{master_updates}**",
        f"- PLACEHOLDERs replaced:      **{placeholder_replaced}**",
        "",
        "## Final master list status (847 codes)",
        "",
        f"- Codes with real paper links:    **{with_real}**",
        f"- Codes still PLACEHOLDER:        **{placeholders}**",
        f"- Codes with no Paper line at all: **{no_paper_line}**",
        "",
    ]
    if unmatched_pdfs:
        rep_lines += ["## Unmatched code directories", ""]
        for cd, files in unmatched_pdfs:
            rep_lines.append(f"- `{cd}`: {len(files)} pdf(s) — _no master entry found_")
        rep_lines.append("")
    REPORT.write_text('\n'.join(rep_lines))
    print(f"Wrote report → {REPORT}")
    print()
    print(f"Final: {with_real}/{total} have real papers; {placeholders} placeholders left; {no_paper_line} entries lack Paper line.")


if __name__ == '__main__':
    main()
