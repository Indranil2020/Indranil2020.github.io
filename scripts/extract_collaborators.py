#!/usr/bin/env python3
"""
Extract collaborators from publications list
Based on your CV and publication data
"""

import json
from collections import Counter
from typing import Dict, List

# Your publications and co-authors (extracted from CV)
from smart_merge_publications import extract_publications_from_html

# Dynamic extraction from publications.html
def get_publications_data():
    try:
        # returns publications, nav_html, profile_links_html, footer_html
        pubs, _, _, _ = extract_publications_from_html('publications.html')
        print(f"Extracted {len(pubs)} publications for collaborator analysis")
        return pubs
    except Exception as e:
        print(f"Error extracting publications: {e}")
        return []

# Collaborator details based on research history and institutions
COLLABORATOR_INFO = {
    "D. P. Samajdar": {
        "name": "Dr. Dip Prakash Samajdar",
        "institution": "PDPM IIITDM Jabalpur",
        "location": "Jabalpur, India",
        "lat": 23.1815, "lon": 79.9864,
        "role": "Ph.D. Supervisor & Primary Collaborator",
        "papers": 0  # Will be counted
    },
    "N. Jain": {
        "name": "Neelesh Jain",
        "institution": "PDPM IIITDM Jabalpur",
        "location": "Jabalpur, India",
        "lat": 23.1815, "lon": 79.9864,
        "role": "Research Collaborator",
        "papers": 0
    },
    "T. D. Das": {
        "name": "Dr. T. D. Das",
        "institution": "NIT Arunachal Pradesh",
        "location": "Yupia, India",
        "lat": 27.1328, "lon": 93.6199,
        "role": "M.Tech Supervisor",
        "papers": 0
    },
    "S. Chakrabarti": {
        "name": "Prof. Subhananda Chakrabarti",
        "institution": "IIT Bombay",
        "location": "Mumbai, India",
        "lat": 19.1334, "lon": 72.9133,
        "role": "Research Collaborator",
        "papers": 0
    },
    "T. Hidouri": {
        "name": "Dr. Tarek Hidouri",
        "institution": "University of Monastir",
        "location": "Monastir, Tunisia",
        "lat": 35.7772, "lon": 10.8261,
        "role": "International Collaborator",
        "papers": 0
    },
    "B. K. Mani": {
        "name": "Prof. Brajesh Kumar Mani",
        "institution": "IIT Delhi",
        "location": "New Delhi, India",
        "lat": 28.5450, "lon": 77.1920,
        "role": "Postdoc Supervisor",
        "papers": 0
    },
    "A. J. Peter": {
        "name": "Dr. A. John Peter",
        "institution": "Madras Christian College",
        "location": "Chennai, India",
        "lat": 13.0338, "lon": 80.2497,
        "role": "Research Collaborator",
        "papers": 0
    },
    "S. Singh": {
        "name": "Sadhna Singh",
        "institution": "PDPM IIITDM Jabalpur",
        "location": "Jabalpur, India",
        "lat": 23.1815, "lon": 79.9864,
        "role": "Research Collaborator",
        "papers": 0
    },
    "D. Roy": {
        "name": "Debamita Roy",
        "institution": "PDPM IIITDM Jabalpur",
        "location": "Jabalpur, India",
        "lat": 23.1815, "lon": 79.9864,
        "role": "Research Collaborator",
        "papers": 0
    },
    "V. Tiwari": {
        "name": "Vikas Tiwari",
        "institution": "PDPM IIITDM Jabalpur",
        "location": "Jabalpur, India",
        "lat": 23.1815, "lon": 79.9864,
        "role": "Research Collaborator",
        "papers": 0
    }
}


def extract_collaborators():
    """Extract and count collaborators"""
    author_counts = Counter()

    # Count papers for each author
    # Name normalization map
    NAME_MAP = {
        "Dip Prakash Samajdar": "D. P. Samajdar",
        "DP Samajdar": "D. P. Samajdar",
        "Dr. Dip Prakash Samajdar": "D. P. Samajdar",
        "Sadhna Singh": "S. Singh",
        "S. Singh": "S. Singh",
        "TD Das": "T. D. Das",
        "T. D. Das": "T. D. Das",
        "N. Jain": "N. Jain",
        "Neelesh Jain": "N. Jain",
        "Subhananda Chakrabarti": "S. Chakrabarti",
        "Prof. Subhananda Chakrabarti": "S. Chakrabarti",
        "Dr. Tarek Hidouri": "T. Hidouri",
        "T. Hidouri": "T. Hidouri",
        "Hidouri Tarek": "T. Hidouri",
        "Mohd Zeeshan": "Mohd Zeeshan", # Not in COLLABORATOR_INFO? Add if causing issues? No, only known ones.
        "B. K. Mani": "B. K. Mani",
        "Prof. Brajesh Kumar Mani": "B. K. Mani",
        "Debamita Roy": "D. Roy",
        "D. Roy": "D. Roy"
    }

    # Count papers for each author
    pubs = get_publications_data()
    for pub in pubs:
        # Authors format in HTML is "Name Name and Name Name" or "Name Name, Name Name"
        # We need to robustly split them.
        author_str = pub.get('authors', '')
        # Remove HTML Highlighting
        author_str = author_str.replace('<span class="me">', '').replace('</span>', '')
        
        # Split by ' and ' or ', '
        if ' and ' in author_str:
            authors = author_str.split(' and ')
        else:
            authors = author_str.split(', ')
            
        for author in authors:
            author = author.strip()
            if not author: continue
            
            # Normalize Name (e.g. "I. Mal" vs "Indranil Mal")
            if "Mal" in author and ("I." in author or "Indranil" in author): 
                continue # Exclude self
            
            # Apply name mapping
            author = NAME_MAP.get(author, author)
                
            author_counts[author] += 1

    # Update paper counts
    collaborators = []
    # Update paper counts strictly for KNOWN collaborators
    collaborators = []
    
    # Reset counts in COLLABORATOR_INFO first to avoid accumulation if script logic changes
    for key in COLLABORATOR_INFO:
        COLLABORATOR_INFO[key]['papers'] = 0

    # Map normalized names to standard keys
    # We need to find which key corresponds to the count
    
    # Ideally, we iterate through our known collaborators and check counts
    # But author_counts has various names.
    
    # Let's iterate author_counts and match to COLLABORATOR_INFO
    for author_name, count in author_counts.items():
        # Check if this author is a known key
        if author_name in COLLABORATOR_INFO:
            COLLABORATOR_INFO[author_name]['papers'] += count
    
    # Now build the list from COLLABORATOR_INFO, excluding 0 papers if desired?
    # User said "change the count nothing else", implying we keep everyone even with 0?
    # Or maybe keep everyone who was there.
    
    for key, info in COLLABORATOR_INFO.items():
        # Only add valid collaborators
        # If papers > 0 or if we want to show them regardless?
        # Let's show them if they have valid paper count usually.
        # But for "strict" check, let's include them.
        
        # Make a copy to return
        c = info.copy()
        collaborators.append(c)

    return sorted(collaborators, key=lambda x: x["papers"], reverse=True)


def main():
    """Main function"""
    collaborators = extract_collaborators()

    data = {
        "total_collaborators": len(collaborators),
        "total_institutions": len(set(c["institution"] for c in collaborators)),
        "collaborators": collaborators
    }

    # Save to JSON
    with open('data/collaborators.json', 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    print(f"Extracted {len(collaborators)} collaborators")
    print(f"Top 5 collaborators:")
    for i, collab in enumerate(collaborators[:5], 1):
        print(f"{i}. {collab['name']}: {collab['papers']} papers")


if __name__ == "__main__":
    main()
