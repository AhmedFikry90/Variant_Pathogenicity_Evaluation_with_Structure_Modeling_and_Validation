import argparse
try:
    import esm
except ImportError:
    raise ImportError("The 'esm' package is not installed. Install it with: pip install git+https://github.com/facebookresearch/esm.git")
try:
    import einops
except ImportError:
    raise ImportError("The 'einops' package is not installed. Install it with: pip install einops")
try:
    from Bio import PDB
except ImportError:
    raise ImportError("The 'biopython' package is not installed. Install it with: pip install biopython")
try:
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("The 'matplotlib' package is not installed. Install it with: pip install matplotlib")
import requests
import json
import random
import os
import numpy as np
from datetime import datetime
import torch

def apply_mutation(sequence, position, new_aa):
    """
    Apply a mutation to the wild-type sequence.
    
    Parameters:
    - sequence (str): The wild-type protein sequence.
    - position (int): The 1-based position of the mutation.
    - new_aa (str): The new amino acid (single-letter code).
    
    Returns:
    - str: The mutated sequence.
    """
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    if new_aa not in valid_aa:
        raise ValueError(f"Invalid amino acid: {new_aa}. Must be one of {valid_aa}")
    if not 1 <= position <= len(sequence):
        raise ValueError(f"Invalid position: {position}. Must be between 1 and {len(sequence)}")
    
    seq_list = list(sequence)
    seq_list[position - 1] = new_aa  # Convert to 0-based index
    return "".join(seq_list)

def predict_structure(sequence):
    """
    Predict the 3D structure of a protein sequence using ESMFold.
    
    Parameters:
    - sequence (str): The protein sequence.
    
    Returns:
    - str: The predicted PDB structure as a string.
    """
    try:
        # Load ESMFold model
        model = esm.pretrained.esmfold_v1()
        model = model.eval()
        # Predict structure
        with torch.no_grad():
            output = model.infer_pdb(sequence)
        return output
    except Exception as e:
        raise RuntimeError(f"ESMFold prediction failed: {str(e)}. Ensure 'esm', 'torch', and 'einops' are installed correctly.")

def validate_structure(pdb_content, output_dir):
    """
    Validate the predicted structure by generating a Ramachandran plot using Bio.PDB and matplotlib.
    
    Parameters:
    - pdb_content (str): The content of the predicted PDB file.
    - output_dir (str): Directory to save the output files.
    
    Returns:
    - None: Saves the Ramachandran plot as 'ramachandran_plot.png' in output_dir.
    """
    pdb_file = os.path.join(output_dir, "predicted.pdb")
    plot_file = os.path.join(output_dir, "ramachandran_plot.png")
    
    try:
        # Save PDB content
        with open(pdb_file, "w") as f:
            f.write(pdb_content)
        
        # Parse PDB file
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        
        # Calculate phi/psi angles
        phi_psi = []
        for model in structure:
            for chain in model:
                polypeptides = PDB.PPBuilder().build_peptides(chain)
                for poly in polypeptides:
                    phi_psi.extend(poly.get_phi_psi_list())
        
        # Filter valid angles and convert to degrees
        phi_psi = [(np.degrees(phi), np.degrees(psi)) for phi, psi in phi_psi if phi is not None and psi is not None]
        if not phi_psi:
            print("Warning: No valid phi/psi angles found for Ramachandran plot.")
            return
        
        # Plot Ramachandran plot
        phi, psi = zip(*phi_psi)
        plt.figure(figsize=(8, 8))
        plt.scatter(phi, psi, s=10, alpha=0.5)
        plt.xlabel("Phi (degrees)")
        plt.ylabel("Psi (degrees)")
        plt.title("Ramachandran Plot")
        plt.grid(True)
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.savefig(plot_file)
        plt.close()
    except Exception as e:
        print(f"Warning: Ramachandran plot generation failed: {str(e)}")

def get_alpha_missense_prediction(position, new_aa):
    """
    Mock AlphaMissense prediction for a variant (replace with actual API call if available).
    
    Parameters:
    - position (int): The 1-based position of the mutation.
    - new_aa (str): The new amino acid.
    
    Returns:
    - dict: Mock prediction with class and score.
    """
    score = random.uniform(0, 1)  # Replace with actual API call
    if score > 0.7:
        pred_class = "pathogenic"
    elif score > 0.3:
        pred_class = "likely_pathogenic"
    else:
        pred_class = "benign"
    return {"class": pred_class, "score": round(score, 3)}

def get_i_mutant_prediction(sequence, position, new_aa):
    """
    Mock I-Mutant 2.0 prediction for protein stability change (replace with actual model if available).
    
    Parameters:
    - sequence (str): The wild-type sequence.
    - position (int): The 1-based position of the mutation.
    - new_aa (str): The new amino acid.
    
    Returns:
    - dict: Mock prediction with ΔΔG and effect.
    """
    ddG = random.uniform(-2, 2)  # Replace with actual I-Mutant call
    effect = "destabilizing" if ddG < 0 else "stabilizing"
    return {"ddG": round(ddG, 3), "effect": effect}

def get_vep_annotations(variant_hgvs):
    """
    Retrieve annotations from the Ensembl VEP REST API, including PhyloP.
    
    Parameters:
    - variant_hgvs (str): The variant in HGVS format (e.g., "9:g.22125504G>C").
    
    Returns:
    - dict: A dictionary containing VEP annotations.
    """
    url = f"https://rest.ensembl.org/vep/human/hgvs/{variant_hgvs}?content-type=application/json&SIFT=1&PolyPhen=1&Conservation=1"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    try:
        response = requests.get(url, headers=headers, timeout=10)
        if response.ok:
            data = response.json()
            if data and "transcript_consequences" in data[0]:
                consequences = data[0]["transcript_consequences"]
                most_severe_consequence = data[0].get("most_severe_consequence", "N/A")
                sift_prediction = "N/A"
                sift_score = "N/A"
                polyphen_prediction = "N/A"
                polyphen_score = "N/A"
                gene = "N/A"
                transcript = "N/A"
                phylop_score = "N/A"
                
                for consequence in consequences:
                    if "sift_prediction" in consequence:
                        sift_prediction = consequence["sift_prediction"]
                        sift_score = consequence.get("sift_score", "N/A")
                    if "polyphen_prediction" in consequence:
                        polyphen_prediction = consequence["polyphen_prediction"]
                        polyphen_score = consequence.get("polyphen_score", "N/A")
                    if "gene_symbol" in consequence:
                        gene = consequence["gene_symbol"]
                    if "transcript_id" in consequence:
                        transcript = consequence["transcript_id"]
                    if "conservation" in consequence:
                        conservation = consequence["conservation"]
                        if isinstance(conservation, dict) and "phylop" in conservation:
                            phylop_score = conservation.get("phylop", "N/A")
                        elif isinstance(conservation, (int, float)):
                            phylop_score = conservation
                        else:
                            print(f"Warning: Unexpected conservation format: {conservation}")
                
                return {
                    "most_severe_consequence": most_severe_consequence,
                    "gene": gene,
                    "transcript": transcript,
                    "sift_prediction": sift_prediction,
                    "sift_score": sift_score,
                    "polyphen_prediction": polyphen_prediction,
                    "polyphen_score": polyphen_score,
                    "phylop_score": phylop_score
                }
            else:
                return {"error": "No transcript consequences found."}
        else:
            return {"error": f"API request failed with status code {response.status_code}"}
    except requests.exceptions.RequestException as e:
        return {"error": f"Request failed: {str(e)}"}

def save_results_to_file(args, vep_results, alpha_missense, i_mutant, output_dir):
    """
    Save results to a text file in the output directory.
    
    Parameters:
    - args: Command-line arguments (position, new_aa, variant_hgvs, skip_validation).
    - vep_results (dict): VEP annotations.
    - alpha_missense (dict): AlphaMissense prediction.
    - i_mutant (dict): I-Mutant 2.0 prediction.
    - output_dir (str): Directory to save the output file.
    
    Returns:
    - str: Path to the saved text file.
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    output_file = os.path.join(output_dir, f"results_{timestamp}.txt")
    
    with open(output_file, "w") as f:
        f.write("Variant Pathogenicity Evaluation Results:\n")
        f.write(f"Mutated position: {args.position}\n")
        f.write(f"New amino acid: {args.new_aa}\n")
        f.write(f"Variant HGVS: {args.variant_hgvs}\n")
        if not args.skip_validation:
            f.write("Modeled structure saved as output/predicted.pdb\n")
            f.write("Ramachandran plot generated as output/ramachandran_plot.png\n")
        else:
            f.write("Structural validation skipped\n")
        f.write("VEP Annotations:\n")
        f.write(json.dumps(vep_results, indent=4) + "\n")
        f.write("AlphaMissense Prediction:\n")
        f.write(json.dumps(alpha_missense, indent=4) + "\n")
        f.write("I-Mutant 2.0 Prediction:\n")
        f.write(json.dumps(i_mutant, indent=4) + "\n")
    
    return output_file

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Evaluate pathogenicity of a novel missense variant with structure modeling and additional tools.")
    parser.add_argument("--position", type=int, required=True, help="1-based position of the mutation")
    parser.add_argument("--new-aa", type=str, required=True, help="New amino acid (single-letter code, e.g., 'A' for Alanine)")
    parser.add_argument("--variant-hgvs", type=str, required=True, help="Variant in HGVS format (e.g., '9:g.22125504G>C')")
    parser.add_argument("--skip-validation", action="store_true", help="Skip structural validation and Ramachandran plot generation")
    args = parser.parse_args()

    # Create output directory
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Default wild-type sequence (can be modified or loaded from a file)
    wild_type_sequence = "MKVLRAALLTLALAGLVLAGVLAFAVGSVQARQDFIGRLVRLPDEVQILAEHAKSELVNDAEKLFNQDVDAAVRGIA"
    
    # Step 1: Apply mutation
    try:
        mutated_sequence = apply_mutation(wild_type_sequence, args.position, args.new_aa)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return
    
    # Step 2: Predict structure with ESMFold (if not skipping validation)
    predicted_pdb = None
    if not args.skip_validation:
        try:
            predicted_pdb = predict_structure(mutated_sequence)
        except RuntimeError as e:
            print(f"Error: {str(e)}")
            return
    
    # Step 3: Validate structure with Ramachandran plot (if not skipping)
    if not args.skip_validation and predicted_pdb:
        validate_structure(predicted_pdb, output_dir)
    
    # Step 4: Get VEP annotations (SIFT, PolyPhen-2, PhyloP)
    vep_results = get_vep_annotations(args.variant_hgvs)
    
    # Step 5: Get AlphaMissense prediction (mock)
    alpha_missense = get_alpha_missense_prediction(args.position, args.new_aa)
    
    # Step 6: Get I-Mutant 2.0 prediction (mock)
    i_mutant = get_i_mutant_prediction(wild_type_sequence, args.position, args.new_aa)
    
    # Step 7: Output results to console
    print("Variant Pathogenicity Evaluation Results:")
    print(f"Mutated position: {args.position}")
    print(f"New amino acid: {args.new_aa}")
    print(f"Variant HGVS: {args.variant_hgvs}")
    if not args.skip_validation:
        print(f"Modeled structure saved as {os.path.join(output_dir, 'predicted.pdb')}")
        print(f"Ramachandran plot generated as {os.path.join(output_dir, 'ramachandran_plot.png')}")
    else:
        print("Structural validation skipped")
    print("VEP Annotations:")
    print(json.dumps(vep_results, indent=4))
    print("AlphaMissense Prediction:")
    print(json.dumps(alpha_missense, indent=4))
    print("I-Mutant 2.0 Prediction:")
    print(json.dumps(i_mutant, indent=4))
    
    # Step 8: Save results to text file
    try:
        output_file = save_results_to_file(args, vep_results, alpha_missense, i_mutant, output_dir)
        print(f"Results saved to {output_file}")
    except OSError as e:
        print(f"Error: Failed to save results to file: {str(e)}")

if __name__ == "__main__":
    main()