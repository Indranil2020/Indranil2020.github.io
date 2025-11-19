#!/usr/bin/env python3
"""
Extract collaborators from publications list
Based on your CV and publication data
"""

import json
from collections import Counter
from typing import Dict, List

# Your publications and co-authors (extracted from CV)
PUBLICATIONS_DATA = {
    "journals": [
        {"authors": ["N. Jain", "I. Mal", "S. Singh", "D. P. Samajdar"]},
        {"authors": ["M. Zeeshan", "I. Mal", "S. Kumawat", "C. K. Vishwakarma", "B. K. Mani"]},
        {"authors": ["N. Jain", "I. Mal", "T. Hidouri", "D. P. Samajdar"]},
        {"authors": ["I. Mal", "D. P. Samajdar"]},
        {"authors": ["N. Jain", "I. Mal", "D. P. Samajdar", "N. Bagga"]},
        {"authors": ["I. Mal", "D. P. Samajdar"]},
        {"authors": ["I. Mal", "R. K. Mahato", "V. Tiwari", "D. P. Samajdar"]},
        {"authors": ["V. Tiwari", "I. Mal", "S. K. Agnihotri", "D. P. Samajdar"]},
        {"authors": ["I. Mal", "D. P. Samajdar"]},
        {"authors": ["T. Hidouri", "S. Nasr", "I. Mal", "D. P. Samajdar", "F. Saidi", "R. Hamila", "H. Maaref"]},
        {"authors": ["I. Mal", "D. P. Panda", "B. Tongbram", "S. Chakrabarti", "D. P. Samajdar"]},
        {"authors": ["T. Hidouri", "M. Biswas", "I. Mal", "S. Nasr", "S. Chakrabarti", "D. P. Samajdar", "F. Saidi"]},
        {"authors": ["D. Roy", "I. Mal", "D. P. Samajdar"]},
        {"authors": ["I. Mal", "D. P. Samajdar"]},
        {"authors": ["T. Hidouri", "I. Mal", "D. P. Samajdar", "F. Saidi", "T. D. Das"]},
        {"authors": ["I. Mal", "J. Jayarubi", "S. Das", "A. S. Sharma", "A. J. Peter", "D. P. Samajdar"]},
        {"authors": ["I. Mal", "D. P. Panda", "B. Tongbram", "D. P. Samajdar", "S. Chakrabarti"]},
        {"authors": ["A. Hazra", "I. Mal", "D. P. Samajdar", "T. D. Das"]},
        {"authors": ["I. Mal", "D. P. Samajdar", "A. J. Peter"]},
        {"authors": ["I. Mal", "D. P. Samajdar", "T. D. Das"]},
    ],
    "book_chapters": [
        {"authors": ["I. Mal", "N. Jain", "D. P. Samajdar"]},
        {"authors": ["N. Jain", "I. Mal", "D. P. Samajdar"]},
        {"authors": ["C. Rajan", "D. P. Samajdar", "I. Mal"]},
        {"authors": ["I. Mal", "D. P. Samajdar"]},
        {"authors": ["I. Mal", "A. Hazra", "D. P. Samajdar", "T. D. Das"]},
        {"authors": ["A. Basu", "A. Saha", "J. Das", "S. Roy", "S. Mitra", "I. Mal", "S. K. Sarkar"]},
    ],
    "conferences": [
        {"authors": ["I. Mal", "S. Singh", "D. P. Samajdar"]},
        {"authors": ["S. Singh", "I. Mal", "D. P. Samajdar", "K. Dutta"]},
        {"authors": ["A. K. Tenwar", "S. Singh", "I. Mal", "D. P. Samajdar"]},
        {"authors": ["S. Saurabh", "S. Singh", "I. Mal", "D. P. Samajdar"]},
        {"authors": ["I. Mal", "D. P. Samajdar"]},
        {"authors": ["D. Roy", "I. Mal", "D. P. Samajdar"]},
    ]
}

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
    for category in PUBLICATIONS_DATA.values():
        for pub in category:
            for author in pub["authors"]:
                if author != "I. Mal":  # Exclude self
                    author_counts[author] += 1

    # Update paper counts
    collaborators = []
    for author, count in author_counts.most_common():
        if author in COLLABORATOR_INFO:
            collab = COLLABORATOR_INFO[author].copy()
            collab["papers"] = count
            collaborators.append(collab)
        else:
            # Unknown collaborator
            collaborators.append({
                "name": author,
                "institution": "Unknown",
                "location": "Unknown",
                "lat": 0, "lon": 0,
                "role": "Collaborator",
                "papers": count
            })

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
