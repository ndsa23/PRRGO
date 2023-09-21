# 🐈 PageRankeRGO (PRRGO)
A Gene Ontology (GO) Tool by Nicholas D'Sa and Luis Solano

## Features
* Produce GO network visualizations
* Find GO terms based on specified keywords
* Rank important GO terms based on the PageRank algorithm
* Retrieve GEO2R differential gene expression data for each GO term

## Installation
1. Install Node.js (comes with npm, a package manager for Node.js): https://nodejs.org/en/download

2. Download the project via Github and unzip if necessary.

3. Open a shell (Command Prompt, PowerShell, Terminal, or other) and change the working directory to the project folder.

4. Check npm packages (these should be pre-installed already, but here they are just in case)
* `cytoscape-dagre@2.5.0`
* `cytoscape@3.23.0`
* `webpack-cli@5.0.1`
* `webpack@5.78.0`

  You can check installed npm packages with `npm list` and download them with `npm install [package-name]`

5. Install Python: https://www.python.org/downloads/
* Python version: Python 3.9.6 (or newer)

6. We recommend installing the following PyPI packages in a virtualenv ([Install](https://virtualenv.pypa.io/en/latest/installation.html), [User Guide](https://virtualenv.pypa.io/en/latest/user_guide.html), [Or using Conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands)). Activate your virtual environment, then run `pip install [package-name]` for each of the below:
* `Django>=4.0`
* `networkx>=3.0`
* `obonet>=1.0.0`
* `scipy>=1.10`

7. Build database:
* `python manage.py makemigrations`
* `python manage.py migrate`

8. Run the server:
* `python manage.py runserver`
* Open http://127.0.0.1:8000/go-viz in a web browser

## Instructions
1. Fill out the form with keyword to search for in GO terms or description, number of GO terms to filter by PageRank, and GEO2R DESeq File(s)
2. Click the "Submit" button. After 15-30 seconds a visualization should appear.
3. The size of the GO term nodes is mapped to PageRank score and the color reflects the GO aspect. Click on a GO term node to view the relevant information for a GO term and the DEG data.
4. (Optional) Download the GO network graph visualization JPEG file or the GO DEG CSV file by clicking the buttons near the bottom of the page. Note that the JPEG file produced will show all nodes (even those offscreen) and will match the current arrangement of nodes rather than the starting arrangement.

To see documentation for key functions, open the `.html` files in the `docs` folder.
