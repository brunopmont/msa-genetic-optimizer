"""
generate_base.py
Data preparation script for generating base alignments or handling raw 
dataset inputs for the MSA optimization pipeline.
"""

import os
import subprocess


def generate_base_alignments() -> None:
    """Generates initial seed alignments for benchmark instances using KAlign."""
    instances = [
        "BB11001", "BB11010", "BB12012", "BB12022", "BB20010", 
        "BB20030", "BB30016", "BB30024", "BB40011", "BB50014"
    ]
    
    os.makedirs("base", exist_ok=True)
    print("[INFO] Starting base alignment generation using KAlign...")
    
    for instance_id in instances:
        input_path = f"bench/bench1.0/bali3/in/{instance_id}"
        output_path = os.path.join("base", f"{instance_id}_base.fasta")
        
        if not os.path.exists(input_path):
            print(f"[WARNING] Input file not found for instance {instance_id}: {input_path}")
            continue
            
        cmd = ["kalign", "-i", input_path, "-o", output_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"[SUCCESS] Base alignment generated for {instance_id}")
        else:
            print(f"[ERROR] Failed to generate base alignment for {instance_id}: {result.stderr.strip()}")


if __name__ == "__main__":
    generate_base_alignments()