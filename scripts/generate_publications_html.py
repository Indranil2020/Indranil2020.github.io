#!/usr/bin/env python3
"""
Generate publications.html from JSON data
"""

import json
from datetime import datetime
from typing import List, Dict


def load_publications(filename: str = 'data/publications.json') -> Dict:
    """Load publications from JSON"""
    with open(filename, 'r', encoding='utf-8') as f:
        return json.load(f)


def categorize_publications(publications: List[Dict]) -> Dict:
    """Categorize publications by type"""
    journals = []
    conferences = []
    book_chapters = []

    for pub in publications:
        venue = pub.get('venue', '').lower()
        title = pub.get('title', '').lower()

        # Categorization logic
        if any(keyword in venue for keyword in ['proceedings', 'conference', 'workshop', 'materials today']):
            conferences.append(pub)
        elif any(keyword in venue for keyword in ['springer proceedings', 'lecture notes', 'advances in']):
            book_chapters.append(pub)
        else:
            journals.append(pub)

    return {
        'journals': sorted(journals, key=lambda x: int(x.get('year', 0)), reverse=True),
        'conferences': sorted(conferences, key=lambda x: int(x.get('year', 0)), reverse=True),
        'book_chapters': sorted(book_chapters, key=lambda x: int(x.get('year', 0)), reverse=True)
    }


def generate_pub_html(pub: Dict, index: int) -> str:
    """Generate HTML for a single publication"""
    authors = pub.get('authors', '').replace('I Mal', '<span class="me">I. Mal</span>')
    title = pub.get('title', '')
    venue = pub.get('venue', '')
    year = pub.get('year', '')
    citations = pub.get('citations', 0)
    doi = pub.get('doi', '')

    html = f'''
            <div class="pub-item">
                <div class="pub-number">[{index}]</div>
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


def generate_html(data: Dict) -> str:
    """Generate complete publications.html"""
    categorized = categorize_publications(data['publications'])

    html_template = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Publications | Indranil Mal</title>
    <!-- Include full stylesheet here -->
</head>
<body>
    <div class="container">
        <h1>Publications</h1>
        <p>Last updated: {last_updated}</p>
        <p>Total Publications: {total} | Total Citations: {citations}</p>

        <section id="journals">
            <h2>Journal Articles ({journal_count})</h2>
            {journals_html}
        </section>

        <section id="book-chapters">
            <h2>Book Chapters ({book_count})</h2>
            {books_html}
        </section>

        <section id="conferences">
            <h2>Conference Publications ({conf_count})</h2>
            {conferences_html}
        </section>
    </div>
</body>
</html>'''

    # Generate HTML for each category
    journals_html = '\\n'.join([generate_pub_html(p, i+1) for i, p in enumerate(categorized['journals'])])
    books_html = '\\n'.join([generate_pub_html(p, i+1) for i, p in enumerate(categorized['book_chapters'])])
    conferences_html = '\\n'.join([generate_pub_html(p, i+1) for i, p in enumerate(categorized['conferences'])])

    total_citations = sum(p.get('citations', 0) for p in data['publications'])

    return html_template.format(
        last_updated=datetime.fromisoformat(data['last_updated']).strftime('%B %d, %Y'),
        total=data['total_publications'],
        citations=total_citations,
        journal_count=len(categorized['journals']),
        book_count=len(categorized['book_chapters']),
        conf_count=len(categorized['conferences']),
        journals_html=journals_html,
        books_html=books_html,
        conferences_html=conferences_html
    )


def main():
    """Main function"""
    print("Generating publications.html from JSON data...")

    # Load data
    data = load_publications('data/publications.json')

    # Generate HTML
    html = generate_html(data)

    # Save
    with open('publications_auto.html', 'w', encoding='utf-8') as f:
        f.write(html)

    print("Generated publications_auto.html")


if __name__ == "__main__":
    main()
