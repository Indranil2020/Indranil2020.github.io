#!/usr/bin/env python3
"""
Intelligent Publication Generator
Properly categorizes and formats publications matching the original structure
"""

import json
import re
from datetime import datetime
from typing import List, Dict, Tuple


class IntelligentPublicationGenerator:
    """Intelligently generates publications.html with proper categorization"""

    # Keywords for intelligent categorization
    CONFERENCE_KEYWORDS = [
        'conference', 'proceedings', 'IEEE', 'workshop', 'symposium',
        'EDKCON', 'international conference', 'national conference'
    ]

    BOOK_KEYWORDS = [
        'springer', 'lecture notes', 'advances in', 'book chapter',
        'handbook', 'encyclopedia', 'proceedings in physics'
    ]

    JOURNAL_KEYWORDS = [
        'journal', 'letters', 'physical review', 'applied physics',
        'materials', 'physics', 'optics', 'electronics', 'nano',
        'science', 'nature', 'ACS', 'RSC', 'Elsevier', 'Wiley',
        'solar energy', 'semiconductor', 'photonics'
    ]

    def __init__(self):
        self.publications = []
        self.categorized = {'journals': [], 'conferences': [], 'book_chapters': []}

    def load_data(self, filename: str = 'data/publications.json') -> Dict:
        """Load publications from JSON"""
        with open(filename, 'r', encoding='utf-8') as f:
            return json.load(f)

    def categorize_publication(self, pub: Dict) -> str:
        """Intelligently categorize publication based on venue"""
        venue = pub.get('venue', '').lower()
        title = pub.get('title', '').lower()

        # Check for arXiv preprints
        if 'arxiv' in venue or 'arxiv' in title:
            return 'preprint'

        # Check for book chapters
        for keyword in self.BOOK_KEYWORDS:
            if keyword in venue:
                return 'book_chapter'

        # Check for conferences
        for keyword in self.CONFERENCE_KEYWORDS:
            if keyword in venue or keyword in title:
                return 'conference'

        # Check for journals
        for keyword in self.JOURNAL_KEYWORDS:
            if keyword in venue:
                return 'journal'

        # Default: if it has a volume/issue pattern, it's likely a journal
        if re.search(r'vol\.?\s*\d+|volume\s*\d+', venue, re.IGNORECASE):
            return 'journal'

        # Default to journal if uncertain
        return 'journal'

    def categorize_all(self, publications: List[Dict]) -> Dict:
        """Categorize all publications"""
        categorized = {
            'journals': [],
            'conferences': [],
            'book_chapters': [],
            'preprints': []
        }

        for pub in publications:
            category = self.categorize_publication(pub)
            if category == 'journal':
                categorized['journals'].append(pub)
            elif category == 'conference':
                categorized['conferences'].append(pub)
            elif category == 'book_chapter':
                categorized['book_chapters'].append(pub)
            elif category == 'preprint':
                categorized['preprints'].append(pub)

        # Sort each category by year (descending)
        for category in categorized:
            categorized[category] = sorted(
                categorized[category],
                key=lambda x: int(x.get('year', 0)),
                reverse=True
            )

        return categorized

    def format_authors(self, authors_str: str) -> str:
        """Format authors with highlighting for I. Mal"""
        # Highlight variations of the author's name
        variations = ['I Mal', 'Indranil Mal', 'I. Mal', 'Mal I']

        for var in variations:
            authors_str = authors_str.replace(var, '<span class="me">I. Mal</span>')

        return authors_str

    def generate_pub_html(self, pub: Dict, number: int) -> str:
        """Generate HTML for single publication with proper formatting"""
        authors = self.format_authors(pub.get('authors', ''))
        title = pub.get('title', '')
        venue = pub.get('venue', '')
        year = pub.get('year', '')
        citations = pub.get('citations', 0)
        doi = pub.get('doi', '')

        html = f'''            <div class="pub-item">
                <div class="pub-number">[{number}]</div>
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

    def read_template(self) -> Tuple[str, str, str]:
        """Read the current publications.html as template"""
        try:
            with open('publications.html', 'r', encoding='utf-8') as f:
                content = f.read()

            # Extract header (everything before journals section)
            header_match = re.search(r'(.*?)<section id="journals">', content, re.DOTALL)
            header = header_match.group(1) if header_match else ''

            # Extract footer (everything after conferences section)
            footer_match = re.search(r'</section>\s*</div>\s*<footer.*', content, re.DOTALL)
            footer = footer_match.group(0) if footer_match else ''

            # Extract stats section
            stats_match = re.search(r'<div class="stats">.*?</div>', content, re.DOTALL)
            stats_section = stats_match.group(0) if stats_match else ''

            return header, stats_section, footer

        except FileNotFoundError:
            return None, None, None

    def generate_complete_html(self, data: Dict) -> str:
        """Generate complete HTML preserving original structure"""
        categorized = self.categorize_all(data['publications'])

        # Read template
        header, stats_section, footer = self.read_template()

        if not header:
            print("Error: Could not read template from publications.html")
            return ""

        # Calculate totals
        total_pubs = len(data['publications'])
        total_citations = sum(p.get('citations', 0) for p in data['publications'])

        # Update stats section
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
                <div class="stat-number">{len(categorized['journals'])}</div>
                <div class="stat-label">Journal Articles</div>
            </div>
            <div class="stat-item">
                <div class="stat-number">{len(categorized['book_chapters'])}</div>
                <div class="stat-label">Book Chapters</div>
            </div>
            <div class="stat-item">
                <div class="stat-number">{len(categorized['conferences'])}</div>
                <div class="stat-label">Conference Papers</div>
            </div>
        </div>'''

        # Generate sections with REVERSE numbering (newest = highest)
        journals_html = self._generate_section_html(categorized['journals'], total_pubs)
        books_html = self._generate_section_html(categorized['book_chapters'], total_pubs - len(categorized['journals']))
        conferences_html = self._generate_section_html(categorized['conferences'],
                                                       total_pubs - len(categorized['journals']) - len(categorized['book_chapters']))

        # Build complete HTML
        html = header
        html += stats_html
        html += f'''

        <!-- Journal Articles -->
        <section id="journals">
            <h2><i class="fas fa-book-open"></i> Journal Articles ({len(categorized['journals'])})</h2>
{journals_html}
        </section>

        <!-- Book Chapters -->
        <section id="book-chapters">
            <h2><i class="fas fa-book"></i> Book Chapters ({len(categorized['book_chapters'])})</h2>
{books_html}
        </section>

        <!-- Conference Publications -->
        <section id="conferences">
            <h2><i class="fas fa-users"></i> Conference Publications ({len(categorized['conferences'])})</h2>
{conferences_html}
        </section>
    </div>

    '''
        html += footer

        return html

    def _generate_section_html(self, pubs: List[Dict], start_number: int) -> str:
        """Generate HTML for a section with reverse numbering"""
        html_parts = []
        current_number = start_number

        for pub in pubs:
            html_parts.append(self.generate_pub_html(pub, current_number))
            current_number -= 1

        return '\n'.join(html_parts)


def main():
    """Main function"""
    print("Intelligent Publication Generator")
    print("=" * 60)

    generator = IntelligentPublicationGenerator()

    # Load data
    print("Loading publications from data/publications.json...")
    data = generator.load_data('data/publications.json')
    print(f"Loaded {data['total_publications']} publications")

    # Generate HTML
    print("Generating publications.html with intelligent categorization...")
    html = generator.generate_complete_html(data)

    if html:
        # Save
        with open('publications.html', 'w', encoding='utf-8') as f:
            f.write(html)

        # Print summary
        categorized = generator.categorize_all(data['publications'])
        print(f"\nSuccessfully generated publications.html!")
        print(f"  - Journals: {len(categorized['journals'])}")
        print(f"  - Book Chapters: {len(categorized['book_chapters'])}")
        print(f"  - Conferences: {len(categorized['conferences'])}")
        print(f"  - Total: {data['total_publications']}")
        print(f"  - Total Citations: {sum(p.get('citations', 0) for p in data['publications'])}")
    else:
        print("ERROR: Failed to generate HTML")


if __name__ == "__main__":
    main()
