# Interactive Tools

This directory contains landing pages for interactive Streamlit-based computational tools.

## Available Tools

### 1. Basis Set Visualizer
- **Page**: `basis-set-viewer.html`
- **Repository**: https://github.com/Indranil2020/DFT_visual
- **Status**: Deployment in progress
- **Description**: Interactive visualization of atomic basis sets from the Basis Set Exchange (BSE)

### 2. k.p Band Structure Calculator
- **Status**: Under development
- **Description**: Calculate and visualize electronic band structures using multi-band k.p perturbation theory

## How to Deploy Streamlit Apps

### Option 1: Streamlit Community Cloud (Recommended - Free)

1. Go to https://streamlit.io/cloud
2. Sign in with GitHub
3. Click "New app"
4. Select repository (e.g., `Indranil2020/DFT_visual`)
5. Specify the main Python file (usually `app.py` or `streamlit_app.py`)
6. Deploy!
7. You'll get a URL like: `https://your-app-name.streamlit.app`
8. Update the link in the corresponding HTML page

### Option 2: Render

1. Go to https://render.com
2. Create a new Web Service
3. Connect your GitHub repository
4. Set build command: `pip install -r requirements.txt`
5. Set start command: `streamlit run app.py --server.port=$PORT --server.address=0.0.0.0`
6. Deploy

### Option 3: Railway

1. Go to https://railway.app
2. Create new project from GitHub repo
3. Add `Procfile` with: `web: streamlit run app.py --server.port=$PORT`
4. Deploy

## Requirements for Each Streamlit App

Each app repository should include:
- `app.py` or `streamlit_app.py` (main application file)
- `requirements.txt` (Python dependencies)
- `README.md` (documentation)
- Optional: `Procfile` for certain deployment platforms

## Updating Links

After deploying a Streamlit app:

1. Get the deployment URL (e.g., `https://basis-set-viewer.streamlit.app`)
2. Edit the corresponding HTML file in this directory
3. Update the "Launch Application" button href
4. Change button class from `btn-secondary` to `btn-primary`
5. Remove the `onclick` alert
6. Commit and push changes

Example:
```html
<!-- Before -->
<a href="#" class="btn btn-secondary" onclick="alert('Coming soon!'); return false;">
    Launch Application (Coming Soon)
</a>

<!-- After -->
<a href="https://your-app.streamlit.app" class="btn btn-primary" target="_blank">
    <i class="fas fa-rocket"></i> Launch Application
</a>
```

## Adding New Tools

To add a new tool:

1. Create a new HTML file in this directory (e.g., `new-tool.html`)
2. Use `basis-set-viewer.html` as a template
3. Update the content with your tool's information
4. Add a card for the tool in `../index.html` under the Tools section
5. Deploy the Streamlit app
6. Update the link in your HTML file

## Contact

For questions or issues: mal@fzu.cz
