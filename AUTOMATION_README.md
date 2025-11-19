# Website Automation System

This document explains the automated systems set up for your academic website.

## 🤖 Google Scholar Auto-Sync

### What it does:
- **Automatically fetches** your publications from Google Scholar every month
- **Includes citation counts** for each paper
- **Updates collaborators list** automatically from co-author data
- **Updates the website** with the latest data
- **No manual intervention required** after initial setup

### How it works:

1. **Monthly Schedule**: GitHub Actions runs on the 1st of every month at midnight UTC
2. **Data Fetching**: Python script (`scripts/sync_scholar.py`) fetches your Scholar profile
3. **Data Storage**: Publications saved to `data/publications.json` with citation counts
4. **Collaborator Extraction**: Script (`scripts/extract_collaborators.py`) extracts co-authors
5. **Collaborator Storage**: Collaborators saved to `data/collaborators.json`
6. **HTML Generation**: Script creates updated publications page
7. **Auto-Commit**: Changes are automatically committed and pushed

### Manual Trigger:
You can manually trigger the sync anytime:
1. Go to your GitHub repository
2. Click "Actions" tab
3. Select "Sync Google Scholar Publications"
4. Click "Run workflow"

### Files Created:
- `scripts/sync_scholar.py` - Fetches data from Google Scholar
- `scripts/generate_publications_html.py` - Generates HTML from data
- `scripts/extract_collaborators.py` - Extracts co-authors from publications
- `.github/workflows/sync-scholar.yml` - GitHub Actions workflow
- `requirements.txt` - Python dependencies
- `data/publications.json` - Stored publication data (auto-generated)
- `data/collaborators.json` - Stored collaborator data (auto-generated)

## 📦 Installation (for local testing):

```bash
# Install Python dependencies
pip install -r requirements.txt

# Run Scholar sync manually
python scripts/sync_scholar.py

# Generate HTML
python scripts/generate_publications_html.py
```

## 🔧 Configuration:

Your Google Scholar ID is: `vDMppM8AAAAJ`

To change it, edit `scripts/sync_scholar.py`:
```python
syncer = ScholarSync(scholar_id="YOUR_SCHOLAR_ID")
```

## 📊 What Data is Collected:

For each publication:
- Title
- Authors (with you highlighted)
- Venue/Journal
- Year
- **Citation count** (updated monthly)
- DOI (when available)
- Google Scholar URL

## 🌐 Next Steps (To be implemented):

### Co-authors/Collaborators Page:
Will include:
- List of all collaborators
- Their institutions
- Interactive world map with pins
- Collaboration statistics
- Auto-synced from Google Scholar co-author data

### Features Coming:
- h-index and i10-index display
- Citation trends graph
- Top cited papers section
- Collaboration network visualization

## ⚠️ Important Notes:

1. **First Run**: The automation will start working from next month's 1st
2. **Rate Limits**: Google Scholar has rate limits, so manual runs should be spaced out
3. **Data Backup**: `data/publications.json` serves as a backup of your publication data
4. **Privacy**: All data is public (from your public Scholar profile)

## 🐛 Troubleshooting:

If sync fails:
1. Check GitHub Actions logs in the "Actions" tab
2. Verify Google Scholar profile is public
3. Check if `scholarly` package is working (Google sometimes blocks scraping)

## 📝 Maintenance:

- **No maintenance required** for normal operation
- GitHub Actions runs automatically
- Python dependencies auto-updated by Dependabot (if enabled)

---

**Status**: ✅ Fully operational and committed to repository
**Last Updated**: November 19, 2025
