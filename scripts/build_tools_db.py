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
import html

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
    """Read full md file text (raw markdown for client-side rendering - legacy)."""
    if not md_rel_path:
        return ''
    p = TOOLS_BASE / md_rel_path
    if not p.exists():
        return ''
    return p.read_text()


# Lazy-init markdown renderer
_md = None
def _get_md():
    global _md
    if _md is None:
        from markdown_it import MarkdownIt
        _md = MarkdownIt('commonmark', {'html': False, 'linkify': True, 'typographer': True})
        _md.enable('table')
        _md.enable('strikethrough')
    return _md


def render_md_to_html(md_rel_path):
    """Render the full md file to HTML at build time (for SEO).
    Returns a tuple (html, plain_text_for_description)."""
    if not md_rel_path:
        return '', ''
    p = TOOLS_BASE / md_rel_path
    if not p.exists():
        return '', ''
    text = p.read_text()
    # Strip the first H1 (we render it in the page header already)
    text_no_h1 = re.sub(r'^#\s+.+?\n', '', text, count=1)
    html = _get_md().render(text_no_h1)
    # Plain-text extraction: strip markdown syntax + html tags for meta description use
    plain = re.sub(r'<[^>]+>', ' ', html)
    plain = re.sub(r'\s+', ' ', plain).strip()
    return html, plain


SECTION_EXTRACT_RE = re.compile(r'^##\s+(.+?)\s*\n+(.+?)(?=\n##\s|\Z)', re.MULTILINE | re.DOTALL)


def extract_sections(md_rel_path):
    """Extract sections from md as {heading_lower: content} for structured data."""
    if not md_rel_path:
        return {}
    p = TOOLS_BASE / md_rel_path
    if not p.exists():
        return {}
    text = p.read_text()
    out = {}
    for m in SECTION_EXTRACT_RE.finditer(text):
        heading = m.group(1).strip().lower()
        body = m.group(2).strip()
        out[heading] = body
    return out


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
        md_rel = e.get('md_link_path', '')
        md_html, md_plain = render_md_to_html(md_rel)
        sections = extract_sections(md_rel)

        # Build paper HTML + JSON-LD papers entries
        papers_jsonld = []
        if e['papers']:
            paper_html = ''.join(
                f'<li><a href="../../scientific_tools_consolidated/{html.escape(p["path"], quote=True)}" target="_blank" rel="noopener"><i class="fas fa-file-pdf"></i> {html.escape(p["name"])}</a></li>'
                for p in e['papers']
            )
            paper_block = f'<h2><i class="fas fa-book"></i> Reference Papers ({len(e["papers"])})</h2><ul class="paper-list">{paper_html}</ul>'
            paper_tag = '<span class="tag"><i class="fas fa-file-pdf"></i> ' + str(len(e['papers'])) + ' paper' + ('s' if len(e['papers']) != 1 else '') + '</span>'
            for p in e['papers']:
                papers_jsonld.append({
                    '@type': 'CreativeWork',
                    'name': p['name'],
                    'url': f"{SITE_URL}/scientific_tools_consolidated/{p['path']}",
                })
        else:
            paper_block = '<h2><i class="fas fa-book"></i> Reference Papers</h2><p class="muted">No paper PDFs uploaded yet for this code.</p>'
            paper_tag = ''

        official = e['official_url'] or ''
        if official.startswith('http'):
            official_block = (
                f'<a href="{html.escape(official, quote=True)}" class="btn-action btn-primary" target="_blank" rel="noopener">'
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

        # Build tagline from overview
        tagline = (e['overview'] or e['note'] or f"{e['name']} - a tool in {cat_clean} / {sub_clean}.")
        tagline = re.sub(r'\s+', ' ', tagline).strip()
        tagline_short = (tagline[:237].rstrip() + '…') if len(tagline) > 240 else tagline

        # Rich meta description: prefer overview + capability summary from md for max content relevance
        capability_text = sections.get('capabilities (critical)', '') or sections.get('capabilities', '') or ''
        capability_text = re.sub(r'[-*]\s*', '', capability_text)
        capability_text = re.sub(r'\s+', ' ', capability_text).strip()
        meta_desc_parts = [tagline]
        if capability_text:
            meta_desc_parts.append('Key capabilities: ' + capability_text[:250])
        meta_desc = ' '.join(meta_desc_parts)
        meta_desc = meta_desc[:297].rstrip() + '…' if len(meta_desc) > 300 else meta_desc

        # Keywords: include tool name variants, section headings, aliases from md
        kw_parts = [
            e['name'], cat_clean, sub_clean,
            e['name'] + ' code', e['name'] + ' software', e['name'] + ' tutorial',
            e['name'] + ' documentation', e['name'] + ' DFT', e['name'] + ' simulation',
            'computational materials science', 'DFT', 'ab initio', 'simulation tool',
            'density functional theory',
        ]
        # Also pull "Scientific domain" bullets from md for richer keywords
        domain_text = ''
        sd_m = re.search(r'\*\*Scientific domain\*\*:\s*([^\n*]+)', md_plain[:2000] if md_plain else '')
        if sd_m:
            domain_text = sd_m.group(1).strip()
            for piece in re.split(r'[,;]', domain_text):
                piece = piece.strip()
                if piece and len(piece) < 60:
                    kw_parts.append(piece)
        keywords = ', '.join(dict.fromkeys([k.strip() for k in kw_parts if k and k.strip()]))

        # Build prev/next within same subcategory
        siblings = entries_by_sub.get((e['category_id'], e['subcategory_id']), [])
        idx_in_sub = next((i for i, x in enumerate(siblings) if x['slug'] == e['slug']), -1)
        prev_e = siblings[idx_in_sub - 1] if idx_in_sub > 0 else None
        next_e = siblings[idx_in_sub + 1] if 0 <= idx_in_sub < len(siblings) - 1 else None
        sub_clean_esc = html.escape(sub_clean, quote=True)
        prev_html = (
            f'<a class="sibling-link prev" href="{html.escape(prev_e["slug"], quote=True)}.html">'
            f'<span class="sib-label"><i class="fas fa-chevron-left"></i> Previous in {sub_clean_esc}</span>'
            f'<span class="sib-name">{html.escape(prev_e["name"])}</span></a>'
        ) if prev_e else '<span class="sibling-link prev disabled"><span class="sib-label">Start of section</span><span class="sib-name">—</span></span>'
        next_html = (
            f'<a class="sibling-link next" href="{html.escape(next_e["slug"], quote=True)}.html">'
            f'<span class="sib-label">Next in {sub_clean_esc} <i class="fas fa-chevron-right"></i></span>'
            f'<span class="sib-name">{html.escape(next_e["name"])}</span></a>'
        ) if next_e else '<span class="sibling-link next disabled"><span class="sib-label">End of section</span><span class="sib-name">—</span></span>'
        siblings_nav = f'<nav class="siblings-nav">{prev_html}{next_html}</nav>'

        # Related tools: list ALL siblings in same subcategory as chips (great internal linking for SEO)
        if siblings:
            related_chips = []
            for sib in siblings:
                cls = 'current' if sib['slug'] == e['slug'] else ''
                if sib['papers']:
                    cls = (cls + ' has-paper').strip()
                chip_cls = f' class="{cls}"' if cls else ''
                related_chips.append(
                    f'            <a href="{html.escape(sib["slug"], quote=True)}.html"{chip_cls}>{html.escape(sib["name"])}</a>'
                )
            related_tools_html = '\n'.join(related_chips)
        else:
            related_tools_html = '            <p class="muted">No sibling tools listed.</p>'

        # Build comprehensive JSON-LD structured data
        canonical_url = f"{SITE_URL}/tools/db/{e['slug']}.html"
        jsonld_main = {
            '@context': 'https://schema.org',
            '@type': 'SoftwareApplication',
            'name': e['name'],
            'alternateName': [e['name'] + ' code', e['name'] + ' software'],
            'applicationCategory': 'ScienceApplication',
            'applicationSubCategory': cat_clean + (' — ' + sub_clean if sub_clean else ''),
            'operatingSystem': 'Linux, macOS, Windows',
            'description': tagline,
            'url': canonical_url,
            'sameAs': [official] if official.startswith('http') else [],
            'author': {'@type': 'Person', 'name': 'Indranil Mal',
                       'url': SITE_URL + '/',
                       'affiliation': {'@type': 'Organization', 'name': 'Institute of Physics, Czech Academy of Sciences'}},
            'keywords': keywords,
            'inLanguage': 'en',
            'isAccessibleForFree': True,
            'mainEntityOfPage': canonical_url,
        }
        if papers_jsonld:
            jsonld_main['citation'] = papers_jsonld

        jsonld_breadcrumb = {
            '@context': 'https://schema.org',
            '@type': 'BreadcrumbList',
            'itemListElement': [
                {'@type': 'ListItem', 'position': 1, 'name': 'Home', 'item': SITE_URL + '/'},
                {'@type': 'ListItem', 'position': 2, 'name': 'Scientific Tools DB', 'item': SITE_URL + '/scientific-tools.html'},
                {'@type': 'ListItem', 'position': 3, 'name': cat_label, 'item': SITE_URL + '/scientific-tools.html#cat=' + str(e['category_id'] or '')},
                {'@type': 'ListItem', 'position': 4, 'name': e['name'], 'item': canonical_url},
            ],
        }
        jsonld_blob = json.dumps(jsonld_main, ensure_ascii=False, indent=2) + '\n</script>\n<script type="application/ld+json">\n' + json.dumps(jsonld_breadcrumb, ensure_ascii=False, indent=2)

        # Escape user-supplied text fields so special chars (& < > ") don't break markup
        def esc(s):
            return html.escape(str(s or ''), quote=True)

        replacements = {
            '{{TITLE}}': esc(f"{e['name']} — {cat_clean} code | Indranil Mal's Tools DB"),
            '{{DESCRIPTION}}': esc(meta_desc),
            '{{KEYWORDS}}': esc(keywords),
            '{{NAME}}': esc(e['name']),
            '{{TAGLINE}}': esc(tagline_short),
            '{{CAT_ID}}': esc(e['category_id'] or ''),
            '{{CAT_LABEL}}': esc(cat_label),
            '{{SUB_LABEL}}': esc(sub_label),
            '{{CAT_COLOR}}': cat_color,  # hex value, safe
            '{{CONFIDENCE}}': esc(e['confidence'] or 'Verified'),
            '{{OVERVIEW}}': esc(e['overview'] or e['note'] or ''),
            '{{OFFICIAL_BLOCK}}': official_block,
            '{{PAPER_BLOCK}}': paper_block,
            '{{PAPER_TAG}}': paper_tag,
            '{{SIBLINGS_NAV}}': siblings_nav,
            '{{RELATED_TOOLS}}': related_tools_html,
            '{{SLUG}}': esc(e['slug']),
            '{{CANONICAL}}': canonical_url,
            '{{MD_HTML}}': md_html or '<p class="muted">No additional documentation file linked yet.</p>',
            '{{JSONLD_BLOB}}': jsonld_blob,
            '{{OG_IMAGE}}': SITE_URL + '/assets/profile.jpg',
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
