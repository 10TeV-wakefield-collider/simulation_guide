# Simulation Guide Documentation

This repository contains the documentation for the 10 TeV Wakefield Collider Design Study simulation tools and examples.

## How to render the documentation on your local computer

### Clone the repository:

   ```bash
   git clone https://github.com/10TeV-wakefield-collider/simulation_guide.git
   cd simulation_guide
   ```

### Setup your conda environment

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

### Build the Documentation Locally

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

## Contributing to the Documentation

### Modifying Documentation

The documentation is written in reStructuredText (RST) format. You can modify the `.rst` files directly in the `docs` directory. For RST syntax reference, see:
- [Sphinx RST Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
- [RST Quick Reference](https://docutils.sourceforge.io/docs/user/rst/quickref.html)

After modifying the RST files, follow the [Building the Documentation Locally](#building-the-documentation-locally) instructions to preview your changes.

### Propose that your changes to be added in the main documentation, via Github

1. Fork the repository, if you have not already done so:
   - Go to [https://github.com/10TeV-wakefield-collider/simulation_guide](https://github.com/10TeV-wakefield-collider/simulation_guide)
   - Click the "Fork" button in the top-right corner

2. Add your fork (replace `YOUR_USERNAME` by your Github username)
   ```bash
   git remote add my_fork https://github.com/YOUR_USERNAME/simulation_guide.git
   cd simulation_guide
   ```

3. Create a new branch for your changes:
   ```bash
   git checkout -b your-branch-name
   ```

4. Make your changes to the documentation files

5. Commit and push your changes:
   ```bash
   git add .
   git commit -m "Description of your changes"
   git push origin your-branch-name
   ```

6. Create a Pull Request:
   - Go to your fork on GitHub
   - Click "New Pull Request"
   - Select your branch
   - Add a description of your changes
   - Submit the pull request

The maintainers will review your changes and provide feedback if needed.
