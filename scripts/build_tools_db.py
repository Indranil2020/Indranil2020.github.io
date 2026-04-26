#!/usr/bin/env python3
"""Build tools_db.json + per-code SEO HTML pages from MASTER_LIST + per-code md files.

This script is the single source of truth for the interactive Scientific Tools page.
Re-run it whenever MASTER_LIST_COMPLETE_372.md or any per-code md file changes.

Usage:
    python3 scripts/build_tools_db.py

Outputs:
    data/tools_db.json                  - Main data file used by the interactive page
    tools/db/<slug>.html                - One SEO page per code (847 of them)
    sitemap.xml                          - Updated with all tool URLs
"""
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MASTER = ROOT / "scientific_tools_consolidated" / "MASTER_LIST_COMPLETE_372.md"
TOOLS_BASE = ROOT / "scientific_tools_consolidated"
OUT_JSON = ROOT / "data" / "tools_db.json"
OUT_HTML_DIR = ROOT / "tools" / "db"
SITEMAP = ROOT / "sitemap.xml"
SITE_URL = "https://indranil2020.github.io"

ENTRY_RE = re.compile(r'^\*\*([0-9]+[a-z]?)\.\s+(.+?)\*\*\*?(?:\s*\(.*?\))?\s*$')
CAT_RE = re.compile(r'^##\s+CATEGORY\s+(\d+):\s+(.+?)\s*\(\d+\s+tools\)\s*$')
SUB_RE = re.compile(r'^###\s+(\d+\.\d+)\s+(.+?)\s*\(\d+\s+tools\)\s*$')
LINK_RE = re.compile(r'^- Link:\s*\[(.+?)\]\((.+?)\)\s*$')
PAPER_RE = re.compile(r'\[([^\]]+\.pdf)\]\(([^)]+\.pdf)\)')
RES_RE = re.compile(r'^- Resources:\s*(.+?)\s*$')
CONF_RE = re.compile(r'^- Confidence:\s*(.+?)\s*$')
NOTE_RE = re.compile(r'^- Notes?:\s*(.+?)\s*$')


def slugify(name):
    s = name.strip()
    s = re.sub(r'[^a-zA-Z0-9._\-]+', '-', s)
    s = re.sub(r'-+', '-', s).strip('-')
    return s or 'unnamed'


def parse_master_list():
    """Walk MASTER_LIST and yield structured entries."""
    text = MASTER.read_text()
    lines = text.split('\n')
    cat_id = cat_name = sub_id = sub_name = None
    entries = []
    i = 0
    while i < len(lines):
        line = lines[i]
        cm = CAT_RE.match(line)
        if cm:
            cat_id, cat_name = cm.group(1), cm.group(2).strip()
            sub_id = sub_name = None
            i += 1
            continue
        sm = SUB_RE.match(line)
        if sm:
            sub_id, sub_name = sm.group(1), sm.group(2).strip()
            i += 1
            continue
        em = ENTRY_RE.match(line)
        if em:
            num, name = em.group(1), em.group(2).strip()
            if 'REMOVED' in name:
                i += 1
                continue
            entry = {
                'num': num,
                'name': name,
                'category_id': cat_id,
                'category': cat_name,
                'subcategory_id': sub_id,
                'subcategory': sub_name,
                'confidence': '',
                'official_url': '',
                'note': '',
                'md_link_text': '',
                'md_link_path': '',
                'papers': [],
            }
            j = i + 1
            while j < len(lines):
                nl = lines[j]
                if ENTRY_RE.match(nl) or CAT_RE.match(nl) or SUB_RE.match(nl):
                    break
                if nl.startswith('---'):
                    break
                m = CONF_RE.match(nl)
                if m:
                    entry['confidence'] = m.group(1)
                m = RES_RE.match(nl)
                if m:
                    entry['official_url'] = m.group(1)
                m = NOTE_RE.match(nl)
                if m:
                    entry['note'] = m.group(1)
                m = LINK_RE.match(nl)
                if m:
                    entry['md_link_text'] = m.group(1)
                    entry['md_link_path'] = m.group(2)
                if nl.lstrip().startswith('- Paper:'):
                    has_placeholder = 'PLACEHOLDER' in nl
                    for pm in PAPER_RE.finditer(nl):
                        entry['papers'].append({
                            'name': pm.group(1),
                            'path': pm.group(2),
                        })
                    entry['paper_placeholder'] = has_placeholder and not entry['papers']
                j += 1
            entries.append(entry)
            i = j
            continue
        i += 1
    return entries


def read_md_overview(md_rel_path):
    """Pull a short overview/abstract from the per-code md file."""
    if not md_rel_path:
        return ''
    p = TOOLS_BASE / md_rel_path
    if not p.exists():
        return ''
    text = p.read_text()
    # Look for ## Overview section, take first paragraph
    m = re.search(r'^##\s+Overview\s*\n+(.+?)(?=\n##\s|\Z)', text, re.MULTILINE | re.DOTALL)
    if m:
        para = m.group(1).split('\n\n')[0].strip()
        para = re.sub(r'\s+', ' ', para)
        return para[:500]
    # Fallback: first non-heading paragraph
    for chunk in text.split('\n\n'):
        c = chunk.strip()
        if c and not c.startswith('#'):
            return re.sub(r'\s+', ' ', c)[:500]
    return ''


def read_md_full_html(md_rel_path):
    """Read full md file text (raw markdown for client-side rendering)."""
    if not md_rel_path:
        return ''
    p = TOOLS_BASE / md_rel_path
    if not p.exists():
        return ''
    return p.read_text()


def build_tree(entries):
    """Hierarchical tree for the mind map: root -> category -> subcategory -> tool."""
    tree = {'name': 'Computational Materials Science', 'children': []}
    cats = {}
    for e in entries:
        cid = e['category_id'] or 'misc'
        cname = e['category'] or 'Misc'
        sid = e['subcategory_id'] or cid
        sname = e['subcategory'] or cname
        c = cats.setdefault(cid, {'id': cid, 'name': cname, 'children': {}})
        s = c['children'].setdefault(sid, {'id': sid, 'name': sname, 'children': []})
        s['children'].append({
            'id': e['slug'],
            'name': e['name'],
            'leaf': True,
            'idx': e['idx'],
        })
    out_cats = []
    for cid in sorted(cats.keys(), key=lambda x: (len(x), x)):
        c = cats[cid]
        subs = []
        for sid in sorted(c['children'].keys(), key=lambda x: (len(x), x)):
            subs.append({
                'id': c['children'][sid]['id'],
                'name': c['children'][sid]['name'],
                'children': c['children'][sid]['children'],
            })
        out_cats.append({'id': c['id'], 'name': c['name'], 'children': subs})
    tree['children'] = out_cats
    return tree


def main():
    print(f"Reading {MASTER}...")
    entries = parse_master_list()
    print(f"Parsed {len(entries)} entries")

    # Add slugs + overview
    for idx, e in enumerate(entries):
        e['slug'] = slugify(e['name'])
        e['idx'] = idx
        e['overview'] = read_md_overview(e.get('md_link_path', ''))

    # Make slugs unique
    seen = {}
    for e in entries:
        s = e['slug']
        if s in seen:
            seen[s] += 1
            e['slug'] = f"{s}-{seen[s]}"
        else:
            seen[s] = 1

    tree = build_tree(entries)

    # Stats
    stats = {
        'total_tools': len(entries),
        'with_papers': sum(1 for e in entries if e['papers']),
        'with_md': sum(1 for e in entries if e['md_link_path']),
        'categories': len(set(e['category_id'] for e in entries if e['category_id'])),
        'subcategories': len(set((e['category_id'], e['subcategory_id']) for e in entries if e['subcategory_id'])),
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        'stats': stats,
        'tree': tree,
        'tools': entries,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {OUT_JSON} ({OUT_JSON.stat().st_size} bytes)")

    # Build per-code SEO pages
    OUT_HTML_DIR.mkdir(parents=True, exist_ok=True)
    template = (ROOT / "scripts" / "tool_page_template.html").read_text()
    written = 0
    for e in entries:
        md_full = read_md_full_html(e.get('md_link_path', ''))
        title = f"{e['name']} | {e['category']} | Computational Tools | Indranil Mal"
        desc = (e['overview'] or e['note'] or f"{e['name']} - {e['category']} / {e['subcategory']} computational materials science tool.")[:300]
        keywords = ', '.join(filter(None, [
            e['name'], e['category'], e['subcategory'],
            'computational materials science', 'DFT', 'simulation tool',
            e['name'] + ' code', e['name'] + ' software', e['name'] + ' tutorial',
        ]))
        # Build paper HTML
        if e['papers']:
            paper_html = ''.join(
                f'<li><a href="../../scientific_tools_consolidated/{p["path"]}" target="_blank" rel="noopener"><i class="fas fa-file-pdf"></i> {p["name"]}</a></li>'
                for p in e['papers']
            )
            paper_block = f'<h2><i class="fas fa-book"></i> Reference Papers</h2><ul class="paper-list">{paper_html}</ul>'
        else:
            paper_block = '<h2><i class="fas fa-book"></i> Reference Papers</h2><p class="muted"><em>No paper PDFs uploaded yet for this code.</em></p>'

        official = e['official_url'] or ''
        official_block = f'<a href="{official}" class="btn-official" target="_blank" rel="noopener"><i class="fas fa-globe"></i> Official Website</a>' if official.startswith('http') else f'<span class="muted">Official link: {official}</span>'

        sub_label = f"{e['subcategory_id']} {e['subcategory']}" if e['subcategory'] else ''
        cat_label = f"{e['category_id']}. {e['category']}" if e['category'] else ''

        # Escape md content for embedding in <script>
        md_escaped = md_full.replace('</script>', '<\\/script>')

        page = template
        page = page.replace('{{TITLE}}', title)
        page = page.replace('{{DESCRIPTION}}', desc)
        page = page.replace('{{KEYWORDS}}', keywords)
        page = page.replace('{{NAME}}', e['name'])
        page = page.replace('{{CATEGORY}}', cat_label)
        page = page.replace('{{SUBCATEGORY}}', sub_label)
        page = page.replace('{{CONFIDENCE}}', e['confidence'])
        page = page.replace('{{OVERVIEW}}', e['overview'] or e['note'] or '')
        page = page.replace('{{OFFICIAL_BLOCK}}', official_block)
        page = page.replace('{{PAPER_BLOCK}}', paper_block)
        page = page.replace('{{SLUG}}', e['slug'])
        page = page.replace('{{CANONICAL}}', f"{SITE_URL}/tools/db/{e['slug']}.html")
        page = page.replace('{{MD_CONTENT}}', md_escaped)

        out_path = OUT_HTML_DIR / f"{e['slug']}.html"
        out_path.write_text(page)
        written += 1
    print(f"Wrote {written} per-code SEO pages to {OUT_HTML_DIR}")

    # Build sitemap
    urls = [
        ('/', '1.0', 'monthly'),
        ('/publications.html', '0.9', 'weekly'),
        ('/scientific-tools.html', '0.95', 'weekly'),
        ('/tools.html', '0.7', 'monthly'),
        ('/collaborators.html', '0.6', 'monthly'),
        ('/resources/fellowships.html', '0.6', 'monthly'),
        ('/tools/basis-set-viewer.html', '0.6', 'monthly'),
    ]
    for e in entries:
        urls.append((f"/tools/db/{e['slug']}.html", '0.7', 'monthly'))

    sm = ['<?xml version="1.0" encoding="UTF-8"?>',
          '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">']
    from datetime import date
    today = date.today().isoformat()
    for path, prio, freq in urls:
        sm.append('  <url>')
        sm.append(f'    <loc>{SITE_URL}{path}</loc>')
        sm.append(f'    <lastmod>{today}</lastmod>')
        sm.append(f'    <changefreq>{freq}</changefreq>')
        sm.append(f'    <priority>{prio}</priority>')
        sm.append('  </url>')
    sm.append('</urlset>')
    SITEMAP.write_text('\n'.join(sm) + '\n')
    print(f"Wrote {SITEMAP} with {len(urls)} URLs")

    print("\n=== Build summary ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")


if __name__ == '__main__':
    main()
