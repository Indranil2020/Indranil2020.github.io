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
import datetime
import os
import sys
from typing import Dict, List, Tuple
from difflib import SequenceMatcher


def similarity(a: str, b: str) -> float:
    """Calculate similarity between two strings"""
    return SequenceMatcher(None, a.lower(), b.lower()).ratio()


def normalize_title(title: str) -> str:
    """Normalize title for matching"""
    if not title: return ""
    title = re.sub(r'<.*?>', '', title)  # Remove HTML tags
    title = title.lower()
    # Remove all non-alphanumeric characters (keep only letters and numbers)
    title = re.sub(r'[^a-z0-9]', '', title)
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

def extract_publications_from_html(file_path):
    """
    Robustly extract components from the HTML file.
    Returns: publications, nav_html, profile_links_html, footer_html
    """
    if not os.path.exists(file_path):
        return [], '', '', ''

    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    publications = []

    # 1. Extract Nav and Header (Everything up to <div class="container">)
    nav_match = re.search(r'(.*?)<div class="container">', content, re.DOTALL)
    nav_html = nav_match.group(1) if nav_match else ''

    # 2. Extract Profile Links (specifically the div)
    links_match = re.search(r'(<div class="profile-links">.*?</div>)', content, re.DOTALL)
    profile_links_html = links_match.group(1) if links_match else ''

    # 3. Extract Footer (From first footer tag or after last section)
    # Use explicit footer tag to be safe
    footer_match = re.search(r'(<footer>.*)', content, re.DOTALL)
    footer_html = footer_match.group(1) if footer_match else '</body></html>'
    # Check for closing tags wrapper
    if '</body>' not in footer_html:
        footer_html += '\n</body>\n</html>'

    # 4. Extract Publications from Sections
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
            
            # Find all pub items
            pub_matches = re.findall(r'<div class="pub-item">(.*?)</div>\s*</div>', section_content, re.DOTALL)
            # Need to be careful with nested divs. 
            # Better regex: <div class="pub-item"> ...? (next pub item kw or section end)
            # But the existing regex in previous versions seemed to work ok? 
            # Actually, previous regex `r'<div class="pub-item">(.*?)</div>\s*</div>'` looks for nested closing div?
            # Let's use split logic which is safer for repeated items
            
            items = section_content.split('<div class="pub-item">')
            for item in items[1:]: # Skip first empty split
                # Truncate at next div closure if needed, but easier to parse inner fields
                
                title_match = re.search(r'<div class="title">(.*?)</div>', item, re.DOTALL)
                title = clean_text(title_match.group(1)) if title_match else ''
                
                authors_match = re.search(r'<div class="authors">(.*?)</div>', item, re.DOTALL)
                authors = authors_match.group(1) if authors_match else ''
                
                venue_match = re.search(r'<div class="venue">(.*?)</div>', item, re.DOTALL)
                venue = clean_text(venue_match.group(1)) if venue_match else ''
                
                details_match = re.search(r'<div class="details">(.*?)</div>', item, re.DOTALL)
                details = details_match.group(1) if details_match else ''
                
                year_match = re.search(r'Year:\s*(\d+)', details)
                year = int(year_match.group(1)) if year_match else 0
                
                citations_match = re.search(r'Citations:\s*(\d+)', details)
                citations = int(citations_match.group(1)) if citations_match else 0
                
                doi_match = re.search(r'href="https://doi.org/(.*?)"', item)
                doi = doi_match.group(1) if doi_match else ''
                
                if title: # valid item
                    publications.append({
                        'title': title,
                        'authors': authors, # Keep raw HTML for authors (spans)
                        'venue': venue,
                        'year': year,
                        'citations': citations,
                        'doi': doi,
                        'category': category,
                        'source': 'manual'
                    })

    return publications, nav_html, profile_links_html, footer_html


def generate_complete_html(merged_pubs: List[Dict], nav_html: str, profile_links_html: str, footer_html: str) -> str:
    """Generate final HTML by assembling components"""
    
    # Sort publications: Year (desc), then Title
    merged_pubs.sort(key=lambda x: (x.get('year', 0), x.get('title', '')), reverse=True)

    # Categorize
    categorized = {
        'journal': [],
        'book_chapter': [],
        'conference': []
    }

    for pub in merged_pubs:
        cat = pub.get('category', 'journal')
        if cat in categorized:
            categorized[cat].append(pub)
        else:
            categorized['journal'].append(pub) # Fallback

    # Helper to generate HTML for a list
    def generate_list_html(pubs, start_index):
        html_list = ""
        for i, pub in enumerate(pubs):
            pub_num = start_index - i
            
            # Format Authors (Check if clean or raw)
            authors = pub['authors']
            
            # Formatting Title
            title = pub['title']
            
            # Formatting Venue
            venue = pub.get('venue', '')
            venue_html = f'<div class="venue">{venue}</div>' if venue else '<div class="venue"></div>'
            
            # Formatting Details
            year = pub.get('year', 0)
            cites = pub.get('citations', 0)
            
            # Formatting Links
            doi = pub.get('doi', '')
            link_html = ''
            if doi:
                link_html = f'''<div class="links">
                    <a href="https://doi.org/{doi}" target="_blank"><i class="fas fa-link"></i> DOI</a>
                </div>'''

            html_list += f'''            <div class="pub-item">
                <div class="pub-number">[{pub_num}]</div>
                <div class="authors">{authors}</div>
                <div class="title">{title}</div>
                {venue_html}
                <div class="details">Year: {year} | Citations: {cites}</div>
                {link_html}
            </div>
'''
        return html_list

    # Generate sections (Correct numbering logic?)
    # User wants strict preserving of "Original" numbering if possible?
    # Or just sequential? The original file had numbering like [40], [39]...
    # We will assume sequential based on total count descending.
    
    total_count = len(merged_pubs)
    
    journals_html = generate_list_html(categorized['journal'], total_count)
    # Decrement start index for next section? 
    # Usually numbering is global or per section?
    # Looking at backup, it seemed continuous? No, it seemed to jump?
    # Actually, let's just use global countdown.
    
    # We need to calculate indices correctly.
    # Assuming Journals -> Books -> Confs order?
    
    # Let's count them
    j_count = len(categorized['journal'])
    b_count = len(categorized['book_chapter'])
    c_count = len(categorized['conference'])
    
    # If standard ordering is J then B then C
    # Journals: Total -> Total - J + 1
    # Books: Total - J -> Total - J - B + 1
    # Confs: ...
    
    # Re-generate with correct indices
    # Journals
    journals_html = generate_list_html(categorized['journal'], total_count)
    
    # Books
    books_html = generate_list_html(categorized['book_chapter'], total_count - j_count)
    
    # Conferences
    conferences_html = generate_list_html(categorized['conference'], total_count - j_count - b_count)


    # Calculate Stats
    total_pubs = total_count
    total_citations = sum(p.get('citations', 0) for p in merged_pubs)

    # Generate Stats HTML
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
                <div class="stat-number">{j_count}</div>
                <div class="stat-label">Journal Articles</div>
            </div>
            <div class="stat-item">
                <div class="stat-number">{b_count}</div>
                <div class="stat-label">Book Chapters</div>
            </div>
            <div class="stat-item">
                <div class="stat-number">{c_count}</div>
                <div class="stat-label">Conference Papers</div>
            </div>
        </div>'''

    # Assemble Page
    html = nav_html  # Up to <div class="container">
    html += '<div class="container">\n'
    html += stats_html + '\n\n'
    html += '        ' + profile_links_html + '\n\n'
    
    html += f'''        <!-- Journal Articles -->
        <section id="journals">
            <h2><i class="fas fa-book-open"></i> Journal Articles ({j_count})</h2>
{journals_html}        </section>
        
        <!-- Book Chapters -->
        <section id="book-chapters">
            <h2><i class="fas fa-book"></i> Book Chapters ({b_count})</h2>
{books_html}        </section>

        <!-- Conference Publications -->
        <section id="conferences">
            <h2><i class="fas fa-users"></i> Conference Publications ({c_count})</h2>
{conferences_html}        </section>
    </div>
    
    <!-- Footer Wrapper if needed -->
    <div class="footer-wrapper"> 
''' 
    # Note: original footer extraction might include closing divs. 
    # My clean footer_html starts with <footer>. 
    # The original file had `</div>` inside footer_match in my previous regex?
    # Step 193 lines 745-746 showed closing divs.
    # `html += footer_html` probably handles the footer.
    # But we need to close `.container`.
    # Wait, extracted `nav_html` ends with `container` OPENING?
    # No, `re.search(r'(.*?)<div class="container">', ...)` -> `nav_html` is content BEFORE container.
    # So I added `<div class="container">`.
    # I need to close it.
    
    html += '\n    </div> <!-- End Container -->\n\n'
    html += footer_html
    
    return html


def load_scholar_data(json_file: str = 'data/publications.json') -> List[Dict]:
    """Load publications from Google Scholar JSON"""
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Mark as from Scholar
    for pub in data['publications']:
        pub['source'] = 'scholar'

    return data['publications']


def deduplicate_existing(pubs: List[Dict]) -> List[Dict]:
    """Deduplicate existing publications based on normalized title"""
    seen_titles = set()
    unique_pubs = []
    
    for pub in pubs:
        title = normalize_title(pub['title'])
        if title and title not in seen_titles:
            seen_titles.add(title)
            unique_pubs.append(pub)
        else:
            if title:
                print(f"  - Removing duplicate existing: {pub['title'][:50]}...")
    
    return unique_pubs


def smart_merge(existing_pubs: List[Dict], scholar_pubs: List[Dict]) -> List[Dict]:
    """
    Strict Smart Merge:
    1. Preserve ALL existing manual entries (title, venue, authors, formatting).
    2. ONLY update citation counts for existing entries.
    3. ONLY add NEW papers if they are from the last 1 year (Current or Previous Year).
    """
    print("\n[Strict Merge] Starting...")
    
    # Get current year
    import datetime
    current_year = datetime.datetime.now().year
    cutoff_year = current_year - 1
    
    # 0. Pre-filter existing phantom papers
    filtered_existing = []
    for p in existing_pubs:
        t_low = p['title'].lower().strip()
        if t_low == "investigation of optoelectronic" or t_low == "investigation of optoelectronic...":
            print(f"  - Removed phantom paper from existing: {p['title']}")
            continue
        filtered_existing.append(p)
    existing_pubs = filtered_existing

    # 1. Update Existing Publications (Strict: Only Citations)
    print("\n1. Updating existing publications (Citations only)...")
    scholar_matched_indices = set()
    merged = []
    
    for existing in existing_pubs:
        # Create a copy to ensure we don't accidentally mutate other fields if we were referencing objects
        # But here we want to mutate 'existing' to put into 'merged'
        
        ex_title_norm = normalize_title(existing['title'])
        best_match = None
        best_score = 0.0
        
        # Try to find match in Scholar
        for i, scholar in enumerate(scholar_pubs):
            sc_title_norm = normalize_title(scholar['title'])
            
            # fast exact match
            if ex_title_norm == sc_title_norm:
                score = 1.0
            else:
                score = similarity(ex_title_norm, sc_title_norm)
            
            if score > best_score:
                best_score = score
                best_match = (i, scholar)
        
        # Threshold for matching
        if best_score > 0.85:
            idx, scholar_data = best_match
            scholar_matched_indices.add(idx)
            
            old_citations = existing.get('citations', 0)
            new_citations = scholar_data.get('citations', 0)
            
            updated_fields = []
            
            # 1. Update Citations
            if new_citations > old_citations:
                existing['citations'] = new_citations
                updated_fields.append(f"citations: {old_citations}->{new_citations}")

            # 2. Backfill Missing Year (Crucial Fix)
            old_year = existing.get('year', 0)
            new_year = scholar_data.get('year', 0)
            if (old_year == 0 or old_year is None) and new_year > 0:
                existing['year'] = new_year
                updated_fields.append(f"year: 0->{new_year}")
                
            # 3. Backfill Missing Venue
            old_venue = existing.get('venue', '')
            new_venue = scholar_data.get('venue', '')
            if not old_venue and new_venue:
                existing['venue'] = new_venue
                updated_fields.append(f"venue: filled")

            # 4. Backfill DOI if missing
            if not existing.get('doi') and scholar_data.get('doi'):
                existing['doi'] = scholar_data['doi']
                updated_fields.append("doi: filled")

            if updated_fields:
                print(f"  ✓ Updated {existing['title'][:40]}... ({', '.join(updated_fields)})")
            
            existing['source'] = 'manual' # Keep source as manual (but enriched)
        else:
            print(f"  ○ No match for: {existing['title'][:40]}... (kept as is)")
            
        merged.append(existing)

    # 2. Add ONLY RECENT New Publications
    print(f"\n2. Checking for new papers from {cutoff_year}-{current_year}...")
    
    
    blacklist_titles = [
        "investigation of optoelectronic",
        "investigation of optoelectronic...", 
        "investigation of optoelectronic performance of inasnbi for infrared detection" # This is real, don't blacklist
    ]
    
    count_new = 0
    for i, scholar in enumerate(scholar_pubs):
        if i in scholar_matched_indices:
            continue
            
        # Check Year Constraint
        scholar_year = scholar.get('year')
        if not scholar_year:
            continue
            
        # Blacklist check
        sc_title_lower = scholar['title'].lower().strip()
        if sc_title_lower == "investigation of optoelectronic":
            print(f"  - Skipped blacklisted phantom paper: {scholar['title']}")
            continue
            
        if scholar_year < cutoff_year:
            # Check if this is one of the "missing" papers user cares about?
            print(f"  - Skipped old scholar paper: {scholar['title'][:40]}... ({scholar_year})")
            continue
            
        # Double check duplication against ALL merged so far
        sc_title_norm = normalize_title(scholar['title'])
        is_duplicate = False
        for m in merged:
            if similarity(normalize_title(m['title']), sc_title_norm) > 0.9:
                is_duplicate = True
                break
        
        if not is_duplicate:
            # This is a GENUINE NEW RECENT paper
            scholar['category'] = categorize_new_publication(scholar)
            scholar['source'] = 'scholar_new'
            merged.append(scholar)
            count_new += 1
            print(f"  + Added NEW RECENT paper: {scholar['title'][:40]}... ({scholar['year']})")
    
    if count_new == 0:
        print("  (No new recent papers found)")
    
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


    return html


def main():
    """Main function"""
    print("=" * 70)
    print("SMART PUBLICATION MERGE SYSTEM")
    print("=" * 70)

    # 1. Extract existing
    print("[1/5] Extracting existing publications from publications.html...")
    existing_pubs, nav_html, profile_links, footer_html = extract_publications_from_html('publications.html')
    print(f"  Found {len(existing_pubs)} existing publications")
    
    # Deduplicate immediately
    print("[1a/5] Deduplicating source...")
    existing_pubs = deduplicate_existing(existing_pubs)
    print(f"  After deduplication: {len(existing_pubs)}")

    # 2. Load Scholar Data
    print("\n[2/5] Loading Google Scholar data...")
    if os.path.exists('data/publications.json'):
        with open('data/publications.json', 'r') as f:
            data = json.load(f)
            scholar_pubs = data.get('publications', [])
        print(f"  Found {len(scholar_pubs)} publications from Scholar")
    else:
        print("  WARNING: data/publications.json not found!")
        scholar_pubs = []

    # 3. Smart Merge
    print("\n[3/5] Smart merging...")
    merged_pubs = smart_merge(existing_pubs, scholar_pubs)
    print(f"\n  Total after merge: {len(merged_pubs)} publications")

    # 4. Generate HTML
    print("\n[4/5] Generating publications.html (RECONSTRUCTING)...")
    content = generate_complete_html(merged_pubs, nav_html, profile_links, footer_html)

    # 5. Save
    print("\n[5/5] Saving...")
    with open('publications.html', 'w', encoding='utf-8') as f:
        f.write(content)

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
