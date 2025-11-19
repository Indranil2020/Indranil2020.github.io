#!/usr/bin/env python3
"""
Smart Publication Merge System
Merges existing manual publications with new Google Scholar data
- Keeps ALL old publications (never loses data)
- Adds NEW publications from Scholar
- Updates citations for matched publications
- NO HARDCODING - everything dynamic
"""

import json
import re
from typing import Dict, List, Tuple
from difflib import SequenceMatcher


def similarity(a: str, b: str) -> float:
    """Calculate similarity between two strings"""
    return SequenceMatcher(None, a.lower(), b.lower()).ratio()


def normalize_title(title: str) -> str:
    """Normalize title for matching"""
    title = re.sub(r'<.*?>', '', title)  # Remove HTML tags
    title = re.sub(r'[^\w\s]', '', title.lower())
    title = re.sub(r'\s+', ' ', title).strip()
    return title


def clean_text(text: str) -> str:
    """Clean HTML tags from text"""
    text = re.sub(r'<span class="me">(.*?)</span>', r'\1', text)
    text = re.sub(r'<.*?>', '', text)
    return text.strip()


def extract_publications_from_html(html_file: str) -> Tuple[List[Dict], str, str]:
    """
    Extract ALL publications from existing HTML
    Returns: (publications, header_template, footer_template)
    """
    with open(html_file, 'r', encoding='utf-8') as f:
        content = f.read()

    publications = []

    # Extract header (everything before first section)
    header_match = re.search(r'(.*?)<section id="journals">', content, re.DOTALL)
    header = header_match.group(1) if header_match else ''

    # Extract footer (everything after last section)
    footer_match = re.search(r'</section>\s*</div>\s*<footer.*', content, re.DOTALL)
    footer = footer_match.group(0) if footer_match else ''

    # Extract from each section
    sections = {
        'journals': 'journal',
        'book-chapters': 'book_chapter',
        'conferences': 'conference'
    }

    for section_id, category in sections.items():
        section_pattern = f'<section id="{section_id}">(.*?)</section>'
        section_match = re.search(section_pattern, content, re.DOTALL)

        if section_match:
            section_content = section_match.group(1)

            # Find all publications in this section
            pub_pattern = r'<div class="pub-item">(.*?)</div>\s*</div>'
            pub_matches = re.findall(pub_pattern, section_content, re.DOTALL)

            for pub_html in pub_matches:
                # Extract all fields
                title_match = re.search(r'<div class="title">(.*?)</div>', pub_html)
                title = clean_text(title_match.group(1)) if title_match else ''

                authors_match = re.search(r'<div class="authors">(.*?)</div>', pub_html)
                authors = authors_match.group(1) if authors_match else ''

                venue_match = re.search(r'<div class="venue">(.*?)</div>', pub_html)
                venue = venue_match.group(1) if venue_match else ''

                details_match = re.search(r'<div class="details">(.*?)</div>', pub_html)
                details = details_match.group(1) if details_match else ''

                year_match = re.search(r'Year:\s*(\d+)', details)
                year = int(year_match.group(1)) if year_match else 0

                citations_match = re.search(r'Citations:\s*(\d+)', details)
                citations = int(citations_match.group(1)) if citations_match else 0

                doi_match = re.search(r'href="https://doi.org/(.*?)"', pub_html)
                doi = doi_match.group(1) if doi_match else ''

                publications.append({
                    'title': title,
                    'authors': authors,
                    'venue': venue,
                    'year': year,
                    'citations': citations,
                    'doi': doi,
                    'category': category,
                    'source': 'manual'  # From manual curation
                })

    return publications, header, footer


def load_scholar_data(json_file: str = 'data/publications.json') -> List[Dict]:
    """Load publications from Google Scholar JSON"""
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Mark as from Scholar
    for pub in data['publications']:
        pub['source'] = 'scholar'

    return data['publications']


def smart_merge(existing_pubs: List[Dict], scholar_pubs: List[Dict]) -> List[Dict]:
    """
    Smart merge: Keep all old + add new + update citations
    """
    merged = []
    scholar_matched = set()

    print("\nMerging publications...")

    # Process existing publications
    for existing in existing_pubs:
        existing_title = normalize_title(existing['title'])
        best_match = None
        best_score = 0.0

        # Try to match with Scholar data
        for i, scholar in enumerate(scholar_pubs):
            if i in scholar_matched:
                continue

            scholar_title = normalize_title(scholar['title'])
            score = similarity(existing_title, scholar_title)

            if score > best_score:
                best_score = score
                best_match = (i, scholar)

        # If good match found (>80% similarity)
        if best_score > 0.8:
            scholar_idx, scholar_data = best_match
            scholar_matched.add(scholar_idx)

            # Update citations, keep everything else from manual data
            existing['citations'] = scholar_data.get('citations', existing['citations'])
            existing['source'] = 'manual_updated'

            # If manual data is missing venue but Scholar has it, use Scholar's
            if not existing['venue'] and scholar_data.get('venue'):
                existing['venue'] = scholar_data['venue']

            # If manual data is missing DOI but Scholar has it, use Scholar's
            if not existing['doi'] and scholar_data.get('doi'):
                existing['doi'] = scholar_data['doi']

            print(f"  ✓ Updated: {existing['title'][:50]}... (citations: {existing['citations']})")
        else:
            # Keep existing as-is
            print(f"  ○ Kept: {existing['title'][:50]}...")

        merged.append(existing)

    # Add NEW publications from Scholar that weren't matched
    for i, scholar in enumerate(scholar_pubs):
        if i not in scholar_matched:
            # Categorize new publication
            scholar['category'] = categorize_new_publication(scholar)
            scholar['source'] = 'scholar_new'
            merged.append(scholar)
            print(f"  + Added NEW: {scholar['title'][:50]}... ({scholar['category']})")

    return merged


def categorize_new_publication(pub: Dict) -> str:
    """Categorize new publication from Scholar"""
    venue = pub.get('venue', '').lower()
    title = pub.get('title', '').lower()

    # Book chapter keywords
    book_keywords = ['springer', 'lecture notes', 'advances in', 'handbook', 'proceedings in physics']
    for kw in book_keywords:
        if kw in venue:
            return 'book_chapter'

    # Conference keywords
    conf_keywords = ['conference', 'proceedings', 'IEEE', 'workshop', 'symposium', 'EDKCON']
    for kw in conf_keywords:
        if kw in venue or kw in title:
            return 'conference'

    # Default: journal
    return 'journal'


def format_authors(authors_str: str) -> str:
    """Format authors with highlighting"""
    if '<span class="me">' in authors_str:
        return authors_str  # Already formatted

    # Highlight variations
    variations = ['I Mal', 'Indranil Mal', 'I. Mal', 'Mal I']
    for var in variations:
        authors_str = authors_str.replace(var, '<span class="me">I. Mal</span>')

    return authors_str


def generate_pub_html(pub: Dict) -> str:
    """Generate HTML for single publication"""
    authors = format_authors(pub.get('authors', ''))
    title = pub.get('title', '')
    venue = pub.get('venue', '')
    year = pub.get('year', '')
    citations = pub.get('citations', 0)
    doi = pub.get('doi', '')

    html = f'''            <div class="pub-item">
                <div class="pub-number">[NUM]</div>
                <div class="authors">{authors}</div>
                <div class="title">{title}</div>
                <div class="venue">{venue}</div>
                <div class="details">Year: {year} | Citations: {citations}</div>'''

    if doi:
        html += f'''
                <div class="links">
                    <a href="https://doi.org/{doi}" target="_blank"><i class="fas fa-link"></i> DOI</a>
                </div>'''

    html += '''
            </div>'''

    return html


def generate_complete_html(merged_pubs: List[Dict], header: str, footer: str) -> str:
    """Generate complete HTML with NO hardcoding"""

    # Categorize all publications
    categorized = {
        'journal': [],
        'book_chapter': [],
        'conference': []
    }

    for pub in merged_pubs:
        category = pub.get('category', 'journal')
        categorized[category].append(pub)

    # Sort each category by year (descending)
    for category in categorized:
        categorized[category] = sorted(
            categorized[category],
            key=lambda x: (int(x.get('year', 0)), x.get('title', '')),
            reverse=True
        )

    # Calculate dynamic stats
    total_pubs = len(merged_pubs)
    total_citations = sum(p.get('citations', 0) for p in merged_pubs)
    journal_count = len(categorized['journal'])
    book_count = len(categorized['book_chapter'])
    conf_count = len(categorized['conference'])

    # Generate stats HTML (DYNAMIC)
    stats_html = f'''        <div class="stats">
            <div class="stat-item">
                <div class="stat-number">{total_pubs}</div>
                <div class="stat-label">Total Publications</div>
            </div>
            <div class="stat-item">
                <div class="stat-number">{total_citations}</div>
                <div class="stat-label">Total Citations</div>
            </div>
            <div class="stat-item">
                <div class="stat-number">{journal_count}</div>
                <div class="stat-label">Journal Articles</div>
            </div>
            <div class="stat-item">
                <div class="stat-number">{book_count}</div>
                <div class="stat-label">Book Chapters</div>
            </div>
            <div class="stat-item">
                <div class="stat-number">{conf_count}</div>
                <div class="stat-label">Conference Papers</div>
            </div>
        </div>'''

    # Generate sections with REVERSE numbering
    def generate_section(pubs: List[Dict], start_num: int) -> str:
        html_parts = []
        num = start_num
        for pub in pubs:
            pub_html = generate_pub_html(pub)
            pub_html = pub_html.replace('[NUM]', str(num))
            html_parts.append(pub_html)
            num -= 1
        return '\n'.join(html_parts)

    journals_html = generate_section(categorized['journal'], total_pubs)
    books_html = generate_section(categorized['book_chapter'],
                                   total_pubs - journal_count)
    conferences_html = generate_section(categorized['conference'],
                                        total_pubs - journal_count - book_count)

    # Build complete HTML
    html = header
    html += stats_html
    html += f'''

        <!-- Journal Articles -->
        <section id="journals">
            <h2><i class="fas fa-book-open"></i> Journal Articles ({journal_count})</h2>
{journals_html}
        </section>

        <!-- Book Chapters -->
        <section id="book-chapters">
            <h2><i class="fas fa-book"></i> Book Chapters ({book_count})</h2>
{books_html}
        </section>

        <!-- Conference Publications -->
        <section id="conferences">
            <h2><i class="fas fa-users"></i> Conference Publications ({conf_count})</h2>
{conferences_html}
        </section>
    </div>

    '''

    # Update footer stats (DYNAMIC)
    footer = re.sub(r'<span class="footer-stat-number">\d+</span>\s*<span class="footer-stat-label"><i class="fas fa-book"></i> Publications</span>',
                   f'<span class="footer-stat-number">{total_pubs}</span>\n                <span class="footer-stat-label"><i class="fas fa-book"></i> Publications</span>',
                   footer)

    html += footer

    return html


def main():
    """Main function"""
    print("=" * 70)
    print("SMART PUBLICATION MERGE SYSTEM")
    print("=" * 70)

    # 1. Extract existing publications
    print("\n[1/5] Extracting existing publications from publications.html...")
    existing_pubs, header, footer = extract_publications_from_html('publications.html')
    print(f"  Found {len(existing_pubs)} existing publications")

    # 2. Load Scholar data
    print("\n[2/5] Loading Google Scholar data...")
    scholar_pubs = load_scholar_data('data/publications.json')
    print(f"  Found {len(scholar_pubs)} publications from Scholar")

    # 3. Smart merge
    print("\n[3/5] Smart merging...")
    merged_pubs = smart_merge(existing_pubs, scholar_pubs)
    print(f"\n  Total after merge: {len(merged_pubs)} publications")

    # 4. Generate HTML
    print("\n[4/5] Generating publications.html (NO hardcoding)...")
    html = generate_complete_html(merged_pubs, header, footer)

    # 5. Save
    print("\n[5/5] Saving...")
    with open('publications.html', 'w', encoding='utf-8') as f:
        f.write(html)

    # Summary
    print("\n" + "=" * 70)
    print("✓ SUCCESS!")
    print("=" * 70)

    # Count by category
    cats = {'journal': 0, 'book_chapter': 0, 'conference': 0}
    for pub in merged_pubs:
        cats[pub.get('category', 'journal')] += 1

    print(f"\nFinal publication count:")
    print(f"  - Journals: {cats['journal']}")
    print(f"  - Book Chapters: {cats['book_chapter']}")
    print(f"  - Conferences: {cats['conference']}")
    print(f"  - TOTAL: {len(merged_pubs)}")
    print(f"  - Total Citations: {sum(p.get('citations', 0) for p in merged_pubs)}")

    # Count by source
    sources = {}
    for pub in merged_pubs:
        source = pub.get('source', 'unknown')
        sources[source] = sources.get(source, 0) + 1

    print(f"\nSources:")
    for source, count in sources.items():
        print(f"  - {source}: {count}")


if __name__ == "__main__":
    main()
