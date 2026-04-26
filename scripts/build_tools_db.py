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

CAT_COLORS = {
    '1':  '#ff6b6b', '2':  '#feca57', '3':  '#48dbfb', '4':  '#1dd1a1',
    '5':  '#ee5a6f', '6':  '#a29bfe', '7':  '#fd9644', '8':  '#54a0ff',
    '9':  '#5f27cd', '10': '#10ac84',
}


def clean_sub_name(name):
    """Strip leading '1.1_' or '1.1 ' prefix and replace underscores with spaces."""
    if not name:
        return ''
    s = re.sub(r'^\d+(?:\.\d+)?[_ ]+', '', str(name))
    s = s.replace('_', ' ').strip()
    return s


def clean_cat_name(name):
    if not name:
        return ''
    return str(name).replace('_', ' ').strip()


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

    # Index entries by (cat_id, sub_id) for prev/next navigation in same subcategory
    entries_by_sub = {}
    for e in entries:
        key = (e['category_id'], e['subcategory_id'])
        entries_by_sub.setdefault(key, []).append(e)

    written = 0
    for e in entries:
        md_full = read_md_full_html(e.get('md_link_path', ''))
        # Build paper HTML
        if e['papers']:
            paper_html = ''.join(
                f'<li><a href="../../scientific_tools_consolidated/{p["path"]}" target="_blank" rel="noopener"><i class="fas fa-file-pdf"></i> {p["name"]}</a></li>'
                for p in e['papers']
            )
            paper_block = f'<h2><i class="fas fa-book"></i> Reference Papers ({len(e["papers"])})</h2><ul class="paper-list">{paper_html}</ul>'
            paper_tag = '<span class="tag"><i class="fas fa-file-pdf"></i> ' + str(len(e['papers'])) + ' paper' + ('s' if len(e['papers']) != 1 else '') + '</span>'
        else:
            paper_block = '<h2><i class="fas fa-book"></i> Reference Papers</h2><p class="muted">No paper PDFs uploaded yet for this code.</p>'
            paper_tag = ''

        official = e['official_url'] or ''
        if official.startswith('http'):
            official_block = (
                f'<a href="{official}" class="btn-action btn-primary" target="_blank" rel="noopener">'
                f'<i class="fas fa-globe"></i> Official Website <i class="fas fa-external-link-alt"></i></a>'
            )
        else:
            official_block = ''

        # Cleaned labels
        cat_clean = clean_cat_name(e['category'])
        sub_clean = clean_sub_name(e['subcategory'])
        cat_label = f"{e['category_id']}. {cat_clean}" if cat_clean else ''
        sub_label = f"{e['subcategory_id']} {sub_clean}" if sub_clean else cat_label
        cat_color = CAT_COLORS.get(str(e['category_id']), '#3498db')

        tagline = (e['overview'] or e['note'] or f"{e['name']} - a tool in {cat_clean} / {sub_clean}.")
        tagline = re.sub(r'\s+', ' ', tagline).strip()
        if len(tagline) > 240:
            tagline = tagline[:237].rstrip() + '…'

        # Build prev/next within same subcategory
        siblings = entries_by_sub.get((e['category_id'], e['subcategory_id']), [])
        idx_in_sub = next((i for i, x in enumerate(siblings) if x['slug'] == e['slug']), -1)
        prev_e = siblings[idx_in_sub - 1] if idx_in_sub > 0 else None
        next_e = siblings[idx_in_sub + 1] if 0 <= idx_in_sub < len(siblings) - 1 else None
        prev_html = (
            f'<a class="sibling-link prev" href="{prev_e["slug"]}.html">'
            f'<span class="sib-label"><i class="fas fa-chevron-left"></i> Previous in {sub_clean}</span>'
            f'<span class="sib-name">{prev_e["name"]}</span></a>'
        ) if prev_e else '<span class="sibling-link prev disabled"><span class="sib-label">Start of section</span><span class="sib-name">—</span></span>'
        next_html = (
            f'<a class="sibling-link next" href="{next_e["slug"]}.html">'
            f'<span class="sib-label">Next in {sub_clean} <i class="fas fa-chevron-right"></i></span>'
            f'<span class="sib-name">{next_e["name"]}</span></a>'
        ) if next_e else '<span class="sibling-link next disabled"><span class="sib-label">End of section</span><span class="sib-name">—</span></span>'
        siblings_nav = f'<nav class="siblings-nav">{prev_html}{next_html}</nav>'

        # Escape md content for embedding in <script>
        md_escaped = md_full.replace('</script>', '<\\/script>')

        # Updated SEO description with cleaned names
        desc_clean = (tagline if tagline else f"{e['name']} - {cat_clean} / {sub_clean} computational materials science tool.")[:300]

        replacements = {
            '{{TITLE}}': f"{e['name']} | {cat_clean} | Computational Tools | Indranil Mal",
            '{{DESCRIPTION}}': desc_clean,
            '{{KEYWORDS}}': ', '.join(filter(None, [
                e['name'], cat_clean, sub_clean,
                'computational materials science', 'DFT', 'simulation tool',
                e['name'] + ' code', e['name'] + ' software', e['name'] + ' tutorial',
            ])),
            '{{NAME}}': e['name'],
            '{{TAGLINE}}': desc_clean,
            '{{CAT_ID}}': str(e['category_id'] or ''),
            '{{CAT_LABEL}}': cat_label,
            '{{SUB_LABEL}}': sub_label,
            '{{CAT_COLOR}}': cat_color,
            '{{CONFIDENCE}}': e['confidence'] or 'Verified',
            '{{OVERVIEW}}': e['overview'] or e['note'] or '',
            '{{OFFICIAL_BLOCK}}': official_block,
            '{{PAPER_BLOCK}}': paper_block,
            '{{PAPER_TAG}}': paper_tag,
            '{{SIBLINGS_NAV}}': siblings_nav,
            '{{SLUG}}': e['slug'],
            '{{CANONICAL}}': f"{SITE_URL}/tools/db/{e['slug']}.html",
            '{{MD_CONTENT}}': md_escaped,
        }
        page = template
        for key, val in replacements.items():
            page = page.replace(key, val)

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
