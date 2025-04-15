# Simulation Guide Documentation

This repository contains the documentation for the 10 TeV Wakefield Collider Design Study simulation tools and examples.

## Setting up the Environment

### Using Conda

1. If you don't have Conda installed, you can download it from the [official website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

2. Create a new conda environment:
   ```bash
   conda create -n simulation-docs python=3.8
   conda activate simulation-docs
   ```

3. Install the required dependencies from requirements.txt:
   ```bash
   pip install -r requirements.txt
   ```

## Building the Documentation Locally

1. Make sure you have activated the conda environment:
   ```bash
   conda activate simulation-docs
   ```

2. Navigate to the docs directory:
   ```bash
   cd docs
   ```

3. Build the HTML documentation:
   ```bash
   make html
   ```

4. View the documentation:
   - The built documentation will be available in `docs/_build/html/`
   - Open `docs/_build/html/index.html` in your web browser
