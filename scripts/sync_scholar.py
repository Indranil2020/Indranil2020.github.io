#!/usr/bin/env python3
"""
Google Scholar Publication Sync Script
Fetches publications from Google Scholar and updates the website
"""

import json
import os
from datetime import datetime
from typing import List, Dict
import re

try:
    from scholarly import scholarly
    SCHOLARLY_AVAILABLE = True
except ImportError:
    SCHOLARLY_AVAILABLE = False
    print("Warning: scholarly package not available. Install with: pip install scholarly")


class ScholarSync:
    """Syncs publications from Google Scholar"""

    def __init__(self, scholar_id: str = "vDMppM8AAAAJ"):
        self.scholar_id = scholar_id
        self.publications = []

    def fetch_publications(self) -> List[Dict]:
        """Fetch all publications from Google Scholar"""
        if not SCHOLARLY_AVAILABLE:
            print("scholarly package not available, using manual data")
            return []

        try:
            # Search for the author
            author = scholarly.search_author_id(self.scholar_id)
            author = scholarly.fill(author, sections=['publications'])

            publications = []
            for pub in author['publications']:
                # Fill in citation details
                pub_filled = scholarly.fill(pub)

                pub_data = {
                    'title': pub_filled['bib'].get('title', ''),
                    'authors': pub_filled['bib'].get('author', ''),
                    'venue': pub_filled['bib'].get('venue', ''),
                    'year': pub_filled['bib'].get('pub_year', ''),
                    'citations': pub_filled.get('num_citations', 0),
                    'url': pub_filled.get('pub_url', ''),
                    'scholar_url': pub_filled.get('author_pub_id', '')
                }

                # Try to extract DOI from URLs or citation
                doi = self.extract_doi(pub_filled)
                if doi:
                    pub_data['doi'] = doi

                publications.append(pub_data)

            self.publications = publications
            return publications

        except Exception as e:
            print(f"Error fetching from Scholar: {e}")
            return []

    def extract_doi(self, publication: Dict) -> str:
        """Extract DOI from publication data"""
        # Check pub_url for DOI
        if 'pub_url' in publication:
            doi_match = re.search(r'10\.\d{4,}/[^\s]+', publication['pub_url'])
            if doi_match:
                return doi_match.group()

        # Check citation text
        if 'bib' in publication and 'citation' in publication['bib']:
            citation = publication['bib']['citation']
            doi_match = re.search(r'10\.\d{4,}/[^\s,]+', citation)
            if doi_match:
                return doi_match.group()

        return ""

    def save_to_json(self, filename: str = 'data/publications.json'):
        """Save publications to JSON file"""
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        data = {
            'last_updated': datetime.now().isoformat(),
            'scholar_id': self.scholar_id,
            'total_publications': len(self.publications),
            'publications': self.publications
        }

        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

        print(f"Saved {len(self.publications)} publications to {filename}")

    def load_from_json(self, filename: str = 'data/publications.json') -> List[Dict]:
        """Load publications from JSON file"""
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self.publications = data.get('publications', [])
                return self.publications
        except FileNotFoundError:
            print(f"File {filename} not found")
            return []


def main():
    """Main function"""
    print("Google Scholar Sync Script")
    print("=" * 50)

    # Initialize syncer
    syncer = ScholarSync(scholar_id="vDMppM8AAAAJ")

    # Fetch publications
    print("Fetching publications from Google Scholar...")
    publications = syncer.fetch_publications()

    if publications:
        print(f"Fetched {len(publications)} publications")

        # Save to JSON
        syncer.save_to_json('data/publications.json')

        # Print summary
        total_citations = sum(p.get('citations', 0) for p in publications)
        print(f"\nTotal citations: {total_citations}")
        print(f"Average citations per paper: {total_citations / len(publications):.1f}")
    else:
        print("No publications fetched or scholarly not available")
        print("Skipping JSON update to prevent data loss.")
        
        # Only create empty file if it really doesn't exist
        if not os.path.exists('data/publications.json'):
            print("Creating empty publications.json as fallback (first run)...")
            os.makedirs('data', exist_ok=True)
            fallback_data = {
                'last_updated': datetime.now().isoformat(),
                'scholar_id': syncer.scholar_id,
                'total_publications': 0,
                'publications': []
            }
            with open('data/publications.json', 'w', encoding='utf-8') as f:
                json.dump(fallback_data, f, indent=2, ensure_ascii=False)
            print("Created empty publications.json file")


if __name__ == "__main__":
    main()
