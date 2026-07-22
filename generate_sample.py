"""
generate_sample.py
Utility script for sampling and filtering sequence datasets or 
preparing subset instances for quick testing and validation.
"""

import os


def generate_sample_dataset() -> None:
    """Creates a small sample subset from the main benchmark repository for testing."""
    sample_dir = "sample_data"
    os.makedirs(sample_dir, exist_ok=True)
    
    source_instance = os.path.join("base", "BB11001_base.fasta")
    target_sample = os.path.join(sample_dir, "sample_BB11001.fasta")
    
    if not os.path.exists(source_instance):
        print(f"[WARNING] Source file not found: {source_instance}. Ensure base alignments are generated first.")
        return
        
    with open(source_instance, "r", encoding="utf-8") as src:
        lines = src.readlines()
        
    # Takes only the first 2 sequences as a smaller sample subset
    with open(target_sample, "w", encoding="utf-8") as tgt:
        tgt.writelines(lines[:4])
        
    print(f"[SUCCESS] Sample dataset created successfully at: {target_sample}")


if __name__ == "__main__":
    print("=== MSA Genetic Optimizer - Sample Generator ===")
    generate_sample_dataset()