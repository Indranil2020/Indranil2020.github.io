#!/usr/bin/env python3
"""
Intelligent Citation Updater
Preserves manually curated publication data, only updates citation counts from Scholar
"""

import json
import re
from typing import Dict, List
from difflib import SequenceMatcher


def similarity(a: str, b: str) -> float:
    """Calculate similarity between two strings"""
    return SequenceMatcher(None, a.lower(), b.lower()).ratio()


def normalize_title(title: str) -> str:
    """Normalize title for matching"""
    # Remove special characters, extra spaces
    title = re.sub(r'[^\w\s]', '', title.lower())
    title = re.sub(r'\s+', ' ', title).strip()
    return title


def extract_publications_from_html(html_file: str) -> List[Dict]:
    """Extract publications from existing HTML file"""
    with open(html_file, 'r', encoding='utf-8') as f:
        content = f.read()

    publications = []

    # Find all publication items
    pub_pattern = r'<div class="pub-item">(.*?)</div>\s*</div>'
    matches = re.findall(pub_pattern, content, re.DOTALL)

    for match in matches:
        # Extract title
        title_match = re.search(r'<div class="title">(.*?)</div>', match)
        title = title_match.group(1) if title_match else ''

        # Extract authors
        authors_match = re.search(r'<div class="authors">(.*?)</div>', match)
        authors = authors_match.group(1) if authors_match else ''

        # Extract venue
        venue_match = re.search(r'<div class="venue">(.*?)</div>', match)
        venue = venue_match.group(1) if venue_match else ''

        # Extract details (year, citations)
        details_match = re.search(r'<div class="details">(.*?)</div>', match)
        details = details_match.group(1) if details_match else ''

        # Extract current citations
        citations_match = re.search(r'Citations:\s*(\d+)', details)
        current_citations = int(citations_match.group(1)) if citations_match else 0

        # Extract year
        year_match = re.search(r'Year:\s*(\d+)', details)
        year = year_match.group(1) if year_match else ''

        # Extract DOI
        doi_match = re.search(r'href="https://doi.org/(.*?)"', match)
        doi = doi_match.group(1) if doi_match else ''

        # Extract number
        num_match = re.search(r'<div class="pub-number">\[(\d+)\]</div>', match)
        number = int(num_match.group(1)) if num_match else 0

        publications.append({
            'title': title,
            'authors': authors,
            'venue': venue,
            'year': year,
            'citations': current_citations,
            'doi': doi,
            'number': number,
            'html': match
        })

    return publications


def match_publications(html_pubs: List[Dict], scholar_pubs: List[Dict]) -> Dict:
    """Match publications from HTML with Scholar data"""
    matches = {}

    for html_pub in html_pubs:
        html_title = normalize_title(html_pub['title'])
        best_match = None
        best_score = 0.0

        for scholar_pub in scholar_pubs:
            scholar_title = normalize_title(scholar_pub['title'])
            score = similarity(html_title, scholar_title)

            if score > best_score:
                best_score = score
                best_match = scholar_pub

        # Match if similarity > 80%
        if best_score > 0.8:
            matches[html_pub['title']] = {
                'scholar_data': best_match,
                'html_data': html_pub,
                'match_score': best_score
            }

    return matches


def update_html_with_citations(html_file: str, matches: Dict, total_citations: int) -> str:
    """Update HTML file with new citation counts"""
    with open(html_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # Update each publication's citations
    for title, match_data in matches.items():
        new_citations = match_data['scholar_data']['citations']
        html_pub = match_data['html_data']

        # Build old and new details strings
        old_details = f'Year: {html_pub["year"]} | Citations: {html_pub["citations"]}'
        new_details = f'Year: {html_pub["year"]} | Citations: {new_citations}'

        # Replace in content (within details div context)
        content = content.replace(old_details, new_details)

    # Update total citations in stats
    stats_pattern = r'(<div class="stat-number">)\d+(</div>\s*<div class="stat-label">Total Citations</div>)'
    content = re.sub(stats_pattern, rf'\g<1>{total_citations}\g<2>', content)

    # Update total publications count (from 35 to 40 if needed)
    total_pubs = len(matches)
    # Note: This might need manual review as we're only counting matched pubs

    return content


def main():
    """Main function"""
    print("Intelligent Citation Updater")
    print("=" * 60)

    # Load Scholar data
    print("Loading Google Scholar data...")
    with open('data/publications.json', 'r', encoding='utf-8') as f:
        scholar_data = json.load(f)

    scholar_pubs = scholar_data['publications']
    print(f"Found {len(scholar_pubs)} publications from Google Scholar")

    # Extract from HTML
    print("Extracting publications from publications.html...")
    html_pubs = extract_publications_from_html('publications.html')
    print(f"Found {len(html_pubs)} publications in HTML")

    # Match publications
    print("Matching publications...")
    matches = match_publications(html_pubs, scholar_pubs)
    print(f"Matched {len(matches)} publications")

    # Calculate new total citations
    total_citations = sum(m['scholar_data']['citations'] for m in matches.values())
    print(f"Total citations: {total_citations}")

    # Update HTML
    print("Updating publications.html with new citation counts...")
    updated_html = update_html_with_citations('publications.html', matches, total_citations)

    # Backup original
    with open('publications.html.backup', 'w', encoding='utf-8') as f:
        with open('publications.html', 'r', encoding='utf-8') as orig:
            f.write(orig.read())

    # Save updated
    with open('publications.html', 'w', encoding='utf-8') as f:
        f.write(updated_html)

    print(f"\n✓ Successfully updated publications.html!")
    print(f"  - Matched: {len(matches)}/{len(html_pubs)} publications")
    print(f"  - Total citations: {total_citations}")
    print(f"  - Backup saved to: publications.html.backup")

    # Report unmatched
    unmatched = len(html_pubs) - len(matches)
    if unmatched > 0:
        print(f"\n⚠ Warning: {unmatched} publications not matched with Scholar data")


if __name__ == "__main__":
    main()
