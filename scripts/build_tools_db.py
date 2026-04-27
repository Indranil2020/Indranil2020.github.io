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
from datetime import date
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

# Match a DOI embedded in a filename, e.g. '..._10.1103_PhysRevB.95.075146.pdf'.
# We treat the FIRST underscore after the registrant ('10.NNNN') as the DOI '/'.
DOI_IN_FILENAME_RE = re.compile(r'(10\.\d{4,9})_([^/]+?)\.pdf$', re.IGNORECASE)


def filename_to_doi_url(filename, fallback_query=None):
    """Convert a paper-PDF filename to an https://doi.org/... URL when a DOI is
    encoded in the filename, else return a Google Scholar search URL using
    `fallback_query` (or None if no fallback can be built)."""
    m = DOI_IN_FILENAME_RE.search(filename)
    if m:
        # Suffix may itself contain underscores that were originally slashes; we
        # only convert the FIRST underscore (registrant separator). Most journal
        # DOIs only have one slash anyway, and doi.org tolerates a fair amount.
        return f'https://doi.org/{m.group(1)}/{m.group(2)}'
    if fallback_query:
        from urllib.parse import quote_plus
        return f'https://scholar.google.com/scholar?q={quote_plus(fallback_query)}'
    return None
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
                        pname = pm.group(1)
                        ppath = pm.group(2)
                        # Build a clean Scholar fallback query: drop .pdf, use
                        # spaces so the search is human-readable.
                        clean_q = re.sub(r'\.pdf$', '', pname, flags=re.IGNORECASE)
                        clean_q = re.sub(r'[_\-]+', ' ', clean_q).strip()
                        entry['papers'].append({
                            'name': pname,
                            'path': ppath,
                            'doi_url': filename_to_doi_url(
                                pname,
                                fallback_query=f"{entry['name']} {clean_q}",
                            ),
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

        # Build paper HTML + JSON-LD papers entries.
        # PDFs are kept local-only; on the public web each paper title links to
        # its DOI page (or Google Scholar as a fallback for non-DOI items).
        papers_jsonld = []
        if e['papers']:
            paper_lis = []
            for p in e['papers']:
                href = p.get('doi_url') or ''
                if href.startswith('https://doi.org/'):
                    icon = 'fas fa-link'
                    badge = '<span class="paper-source">DOI</span>'
                elif href.startswith('https://scholar.google.com/'):
                    icon = 'fas fa-search'
                    badge = '<span class="paper-source">Scholar</span>'
                else:
                    icon = 'fas fa-book'
                    badge = ''
                if href:
                    paper_lis.append(
                        f'<li><a href="{html.escape(href, quote=True)}" '
                        f'target="_blank" rel="noopener">'
                        f'<i class="{icon}"></i> {html.escape(p["name"])}</a> {badge}</li>'
                    )
                else:
                    paper_lis.append(f'<li>{html.escape(p["name"])}</li>')
            paper_html = ''.join(paper_lis)
            paper_block = (
                f'<h2><i class="fas fa-book"></i> Reference Papers ({len(e["papers"])})</h2>'
                f'<ul class="paper-list">{paper_html}</ul>'
            )
            paper_tag = (
                '<span class="tag"><i class="fas fa-book"></i> '
                + str(len(e['papers'])) + ' paper'
                + ('s' if len(e['papers']) != 1 else '') + '</span>'
            )
            for p in e['papers']:
                jsonld_entry = {
                    '@type': 'ScholarlyArticle',
                    'name': p['name'],
                }
                if p.get('doi_url'):
                    jsonld_entry['url'] = p['doi_url']
                    if p['doi_url'].startswith('https://doi.org/'):
                        jsonld_entry['identifier'] = {
                            '@type': 'PropertyValue',
                            'propertyID': 'DOI',
                            'value': p['doi_url'].replace('https://doi.org/', ''),
                        }
                papers_jsonld.append(jsonld_entry)
        else:
            paper_block = (
                '<h2><i class="fas fa-book"></i> Reference Papers</h2>'
                '<p class="muted">Reference papers are not yet linked for this code.</p>'
            )
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

    # ------------------------------------------------------------------
    # Static directory hub page: /tools/index.html
    # ------------------------------------------------------------------
    # Why we need this in addition to scientific-tools.html (which is the
    # interactive D3 / 3D-graph view): scientific-tools.html populates its
    # tool list from tools_db.json *via JavaScript*, so non-JS crawlers
    # (and Google's mobile crawler with limited JS budget) can't reliably
    # follow links into the 847 per-tool SEO pages from there.
    #
    # This static index lists every code as a plain <a href> link grouped
    # by category and subcategory, giving every per-tool page a strong
    # internal inbound link from one high-authority hub. Also acts as a
    # graceful fallback for users without JavaScript.
    cats_seen = []
    cats = {}
    for e in entries:
        cid = str(e['category_id'] or '')
        cat_clean = clean_cat_name(e['category'])
        sid = str(e['subcategory_id'] or '')
        sub_clean = clean_sub_name(e['subcategory'])
        cat_key = cid
        if cat_key not in cats:
            cats[cat_key] = {
                'id': cid,
                'label': f"{cid}. {cat_clean}" if cat_clean else cid,
                'color': CAT_COLORS.get(cid, '#3498db'),
                'subs': {},
                '_order': [],
            }
            cats_seen.append(cat_key)
        bucket = cats[cat_key]
        if sid not in bucket['subs']:
            bucket['subs'][sid] = {
                'id': sid,
                'label': f"{sid} {sub_clean}".strip() if sub_clean else (cat_clean or 'Other'),
                'tools': [],
            }
            bucket['_order'].append(sid)
        bucket['subs'][sid]['tools'].append(e)

    rows_html = []
    for ckey in cats_seen:
        c = cats[ckey]
        rows_html.append(
            f'<section class="cat-block" id="cat-{html.escape(c["id"], quote=True)}">'
            f'  <h2 style="border-left:6px solid {c["color"]};padding-left:0.6rem;margin-top:1.5rem;">'
            f'{html.escape(c["label"])} '
            f'<small style="color:#666;font-weight:normal;">'
            f'({sum(len(s["tools"]) for s in c["subs"].values())} tools)</small></h2>'
        )
        for skey in c['_order']:
            s = c['subs'][skey]
            tools_sorted = sorted(s['tools'], key=lambda t: t['name'].lower())
            tool_links = ' &middot; '.join(
                f'<a href="db/{html.escape(t["slug"], quote=True)}.html">{html.escape(t["name"])}</a>'
                for t in tools_sorted
            )
            rows_html.append(
                f'  <div class="sub-block"><h3 style="font-size:1.0rem;color:#444;margin-top:0.8rem;margin-bottom:0.3rem;">'
                f'{html.escape(s["label"])} <small style="color:#888;font-weight:normal;">({len(s["tools"])})</small></h3>'
                f'  <p style="line-height:1.9;">{tool_links}</p></div>'
            )
        rows_html.append('</section>')

    total_tools = len(entries)
    cat_summary = ' &middot; '.join(
        f'<a href="#cat-{html.escape(c["id"], quote=True)}">{html.escape(c["label"])}</a>'
        for c in cats.values() if c['id']
    )
    today = date.today().isoformat()
    hub_canonical = f'{SITE_URL}/tools/'
    hub_title = f"Computational Materials Science Codes — Directory of {total_tools} Tools"
    hub_desc = (
        f"Curated directory of {total_tools} computational materials science codes "
        "across DFT, TDDFT, DMFT, GW/BSE, tight-binding, phonons, molecular dynamics, "
        "structure prediction, post-processing and machine-learning frameworks. "
        "Browse by category; each entry links to detailed documentation, official site, "
        "and DOI references."
    )
    hub_jsonld = json.dumps({
        '@context': 'https://schema.org',
        '@type': 'CollectionPage',
        'name': hub_title,
        'description': hub_desc,
        'url': hub_canonical,
        'isPartOf': {'@type': 'WebSite', 'name': "Indranil Mal", 'url': SITE_URL + '/'},
        'mainEntity': {
            '@type': 'ItemList',
            'numberOfItems': total_tools,
            'itemListElement': [
                {
                    '@type': 'ListItem',
                    'position': i + 1,
                    'url': f"{SITE_URL}/tools/db/{e['slug']}.html",
                    'name': e['name'],
                }
                for i, e in enumerate(entries)
            ],
        },
    }, ensure_ascii=False)

    hub_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>{html.escape(hub_title)}</title>
<meta name="description" content="{html.escape(hub_desc, quote=True)}">
<link rel="canonical" href="{hub_canonical}">
<meta property="og:type" content="website">
<meta property="og:title" content="{html.escape(hub_title)}">
<meta property="og:description" content="{html.escape(hub_desc, quote=True)}">
<meta property="og:url" content="{hub_canonical}">
<meta name="twitter:card" content="summary">
<meta name="twitter:title" content="{html.escape(hub_title)}">
<meta name="twitter:description" content="{html.escape(hub_desc, quote=True)}">
<script type="application/ld+json">{hub_jsonld}</script>
<link rel="preconnect" href="https://cdnjs.cloudflare.com">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
<style>
:root {{ --bg:#fff; --fg:#1a2238; --muted:#5a6478; --primary:#2c4ed4; --border:#e2e8f0; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; background:var(--bg); color:var(--fg); line-height:1.55; }}
header.site-nav {{ background:#141c37; color:#fff; padding:0.8rem 1.4rem; display:flex; flex-wrap:wrap; align-items:center; gap:1.2rem; box-shadow:0 2px 6px rgba(0,0,0,0.15); }}
header.site-nav .brand {{ font-weight:700; font-size:1.05rem; color:#fff; text-decoration:none; }}
header.site-nav nav a {{ color:#cbd5e1; text-decoration:none; font-size:0.92rem; margin-right:0.9rem; }}
header.site-nav nav a:hover {{ color:#fff; }}
.wrap {{ max-width:1180px; margin:0 auto; padding:1.4rem 1.4rem 4rem; }}
h1 {{ font-size:1.85rem; margin:0.6rem 0 0.4rem; }}
.lede {{ color:var(--muted); margin:0 0 1rem; max-width:78ch; }}
.cat-jump {{ background:#f4f6fb; border:1px solid var(--border); border-radius:8px; padding:0.7rem 0.9rem; margin:1rem 0 0.5rem; font-size:0.92rem; }}
.cat-jump a {{ color:var(--primary); text-decoration:none; }}
.cat-jump a:hover {{ text-decoration:underline; }}
.cat-block a {{ color:var(--primary); text-decoration:none; }}
.cat-block a:hover {{ text-decoration:underline; }}
.cat-block h2 {{ font-size:1.25rem; }}
.cat-block h3 {{ font-size:1.0rem; }}
.toolbar {{ display:flex; gap:0.8rem; flex-wrap:wrap; margin-bottom:0.8rem; }}
.toolbar input {{ flex:1; min-width:240px; padding:0.55rem 0.8rem; border:1px solid var(--border); border-radius:6px; font:inherit; }}
.toolbar a.btn {{ background:var(--primary); color:#fff; padding:0.55rem 0.9rem; border-radius:6px; text-decoration:none; font-size:0.9rem; }}
.muted {{ color:var(--muted); font-size:0.88rem; }}
footer {{ text-align:center; padding:1.5rem 1rem; color:var(--muted); border-top:1px solid var(--border); font-size:0.85rem; }}
</style>
</head>
<body>
<header class="site-nav">
  <a class="brand" href="/">Indranil Mal</a>
  <nav>
    <a href="/"><i class="fas fa-home"></i> Home</a>
    <a href="/publications.html"><i class="fas fa-file-alt"></i> Publications</a>
    <a href="/scientific-tools.html"><i class="fas fa-project-diagram"></i> Interactive Tools</a>
    <a href="/tools/" aria-current="page"><i class="fas fa-list"></i> Tools Directory</a>
    <a href="/collaborators.html"><i class="fas fa-users"></i> Collaborators</a>
  </nav>
</header>
<div class="wrap">
  <h1>{html.escape(hub_title)}</h1>
  <p class="lede">{html.escape(hub_desc)}</p>
  <div class="toolbar">
    <input id="filter" type="search" placeholder="Filter by code name (e.g. VASP, Quantum ESPRESSO, Z2Pack…)" autocomplete="off">
    <a class="btn" href="/scientific-tools.html"><i class="fas fa-project-diagram"></i> Interactive view</a>
  </div>
  <div class="cat-jump"><strong>Jump to category:</strong> {cat_summary}</div>
  {''.join(rows_html)}
</div>
<footer>
  <p>Curated by <a href="/">Indranil Mal</a> &middot; updated {today} &middot;
     <a href="https://github.com/Indranil2020/Indranil2020.github.io">source on GitHub</a></p>
</footer>
<script>
// Lightweight client-side filter (progressive enhancement; the page is
// already fully crawlable without JS).
(function() {{
  const inp = document.getElementById('filter');
  if (!inp) return;
  const cats = document.querySelectorAll('.cat-block');
  inp.addEventListener('input', () => {{
    const q = inp.value.trim().toLowerCase();
    cats.forEach(cat => {{
      let any = false;
      cat.querySelectorAll('.sub-block').forEach(sb => {{
        let subAny = false;
        sb.querySelectorAll('a').forEach(a => {{
          const hit = !q || a.textContent.toLowerCase().includes(q);
          a.style.display = hit ? '' : 'none';
          // keep the dot separator visually consistent: hide the next text node too
          if (a.nextSibling && a.nextSibling.nodeType === 3) {{
            a.nextSibling.nodeValue = hit ? a.nextSibling.nodeValue : '';
          }}
          if (hit) subAny = true;
        }});
        sb.style.display = subAny ? '' : 'none';
        if (subAny) any = true;
      }});
      cat.style.display = any ? '' : 'none';
    }});
  }});
}})();
</script>
</body>
</html>
"""
    hub_path = OUT_HTML_DIR.parent / 'index.html'
    hub_path.write_text(hub_html)
    print(f"Wrote tools directory hub to {hub_path}")

    # Build sitemap
    urls = [
        ('/', '1.0', 'monthly'),
        ('/publications.html', '0.9', 'weekly'),
        ('/scientific-tools.html', '0.95', 'weekly'),
        ('/tools/', '0.95', 'weekly'),
        ('/tools.html', '0.7', 'monthly'),
        ('/collaborators.html', '0.6', 'monthly'),
        ('/resources/fellowships.html', '0.6', 'monthly'),
        ('/tools/basis-set-viewer.html', '0.6', 'monthly'),
    ]
    for e in entries:
        urls.append((f"/tools/db/{e['slug']}.html", '0.7', 'monthly'))

    sm = ['<?xml version="1.0" encoding="UTF-8"?>',
          '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">']
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
