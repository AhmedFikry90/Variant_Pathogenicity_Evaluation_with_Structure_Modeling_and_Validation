# Variant Pathogenicity Evaluation with Structure Modeling and Validation

This repository contains a Python script for evaluating the pathogenicity of novel missense variants by modeling their protein structures, validating them, and using multiple prediction tools: Ensembl VEP (SIFT, PolyPhen-2, PhyloP), AlphaMissense, and I-Mutant 2.0. Results are saved to a text file in an `output` directory.

## Prerequisites
- Python 3.7+
- Required libraries: `esm`, `biopython`, `requests`, `torch`, `matplotlib`, `einops`

## Installation
1. **Set up a virtual environment** (recommended):
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
2. **Install PyTorch** (CPU or GPU version, see https://pytorch.org/get-started/locally/):
   ```bash
   pip install torch
   ```
   For GPU support (recommended for ESMFold):
   ```bash
   pip install torch --index-url https://download.pytorch.org/whl/cu118
   ```
3. **Install the `esm` package for ESMFold**:
   ```bash
   pip install git+https://github.com/facebookresearch/esm.git
   ```
4. **Install other dependencies**:
   ```bash
   pip install biopython requests matplotlib einops
   ```
5. **Verify installation**:
   ```bash
   python3 -c "import esm; import Bio; import requests; import torch; import matplotlib; import einops; print('All dependencies installed')"
   ```

## Usage
1. Save the script as `variant_pathogenicity_with_structure.py`.
2. Run the script with command-line arguments for mutation position, new amino acid, and HGVS variant:
   ```bash
   python3 variant_pathogenicity_with_structure.py --position 50 --new-aa A --variant-hgvs "9:g.22125504G>C"
   ```
   To skip structural validation:
   ```bash
   python3 variant_pathogenicity_with_structure.py --position 50 --new-aa A --variant-hgvs "9:g.22125504G>C" --skip-validation
   ```
3. Review the outputs:
   - `output/predicted.pdb`: The modeled protein structure (if validation not skipped).
   - `output/ramachandran_plot.png`: The Ramachandran plot for validation (if validation not skipped).
   - `output/results_YYYYMMDD_HHMM.txt`: Text file with all results (VEP, AlphaMissense, I-Mutant 2.0).
   - Annotations printed to the console.

## Command-Line Arguments
- `--position`: The 1-based position of the mutation (e.g., `50`).
- `--new-aa`: The new amino acid (single-letter code, e.g., `A` for Alanine).
- `--variant-hgvs`: The variant in HGVS format (e.g., `9:g.22125504G>C`).
- `--skip-validation`: Skip structural validation and Ramachandran plot generation (optional).

## Output
- **output/predicted.pdb**: The 3D structure of the mutated protein (if validation not skipped).
- **output/ramachandran_plot.png**: A plot showing phi-psi angles for structural validation (if validation not skipped).
- **output/results_YYYYMMDD_HHMM.txt**: A text file containing:
  - Mutation details (position, new amino acid, HGVS variant).
  - Paths to PDB and Ramachandran plot files (if validation not skipped).
  - VEP annotations (gene, transcript, SIFT, PolyPhen-2, PhyloP).
  - AlphaMissense prediction (class and score, mock implementation).
  - I-Mutant 2.0 prediction (ΔΔG and effect, mock implementation).
- **Console Output**: Same information as the text file.

## Troubleshooting
- **ModuleNotFoundError**:
  - For `esm`:
    ```bash
    pip install git+https://github.com/facebookresearch/esm.git
    ```
  - For `einops`:
    ```bash
    pip install einops
    ```
  - For `biopython`, `requests`, or `matplotlib`:
    ```bash
    pip install biopython requests matplotlib
    ```
  - Verify the Python path:
    ```bash
    which python3
    python3 -c "import sys; print(sys.executable)"
    ```
  - Use a virtual environment to avoid conflicts.
- **PyTorch Issues**: Ensure PyTorch is compatible with your system. Check GPU availability:
  ```bash
  python3 -c "import torch; print(torch.cuda.is_available())"
  ```
- **VEP API Errors**: If the Ensembl VEP API returns unexpected data (e.g., `conservation` as a float):
  - The script handles this by assigning float values directly to `phylop_score` or setting to "N/A" if missing.
  - Verify the HGVS format using tools like VariantValidator (https://variantvalidator.org/).
- **Resource Constraints**: ESMFold is memory-intensive. Use a GPU or run in Google Colab:
  - Colab notebook: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ESMFold.ipynb
- **Permission Errors**: Use `--user` or a virtual environment:
  ```bash
  pip install --user biopython requests matplotlib einops
  ```

## Notes
- Ensure an internet connection for the Ensembl VEP REST API.
- Novel variants may lack predictions from some tools (e.g., SIFT, PolyPhen-2, PhyloP). The script handles variable API responses (e.g., `conservation` as a float).
- AlphaMissense and I-Mutant 2.0 use mock implementations; replace with actual APIs or models if available.
- ESMFold requires significant computational resources; consider a GPU or cloud-based environment like Colab.
- The wild-type sequence is hardcoded but can be modified or loaded from a file.
- The `output` directory is created automatically if it doesn't exist.

## License
MIT License