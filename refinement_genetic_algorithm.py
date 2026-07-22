"""
refinement_genetic_algorithm.py
Core implementation of the Genetic Algorithm (GA) with KAlign Warm Start 
and conservative crossover operators for Multiple Sequence Alignment (MSA).
"""

import os
import subprocess
import random
import copy
import re
import time


def read_aligned_fasta(file_path: str) -> tuple[list[str], list[str]]:
    """Reads a FASTA file and returns headers and aligned sequences."""
    headers, sequences, current_seq = [], [], ""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    headers.append(line)
                    if current_seq:
                        sequences.append(current_seq)
                        current_seq = ""
                elif line:
                    current_seq += line
            if current_seq:
                sequences.append(current_seq)
        return headers, sequences
    except FileNotFoundError:
        return [], []


def save_fasta_result(headers: list[str], alignment: list[str], output_path: str) -> None:
    """Saves the final alignment matrix back into FASTA format."""
    with open(output_path, 'w', encoding='utf-8') as f:
        for i in range(len(alignment)):
            f.write(f"{headers[i]}\n{alignment[i]}\n")


def calculate_wsp_fitness(alignment: list[str]) -> float:
    """Calculates the Weighted Sum of Pairs (WSP) score for an alignment matrix."""
    score = 0.0
    length = len(alignment[0])
    num_seqs = len(alignment)
    
    for col in range(length):
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                l1, l2 = alignment[i][col], alignment[j][col]
                if l1 != '-' and l2 != '-' and l1 == l2:
                    score += 1.0
                elif l1 != '-' or l2 != '-':
                    score -= 0.5
    return score


def conservative_crossover(parent_1: list[str], parent_2: list[str]) -> tuple[list[str], list[str]]:
    """Performs topological biological crossover preserving building blocks."""
    length = min(len(parent_1[0]), len(parent_2[0]))
    if length < 3:
        return parent_1, parent_2
        
    cut_point = random.randint(1, length - 1)
    child_1, child_2 = [], []
    
    for i in range(len(parent_1)):
        left_1 = parent_1[i][:cut_point]
        letters_left_1 = len(left_1.replace('-', ''))
        seen, cut_p2 = 0, 0
        
        if letters_left_1 > 0:
            for j, char in enumerate(parent_2[i]):
                if char != '-':
                    seen += 1
                if seen == letters_left_1:
                    cut_p2 = j + 1
                    break
        child_1.append(left_1 + parent_2[i][cut_p2:])
        
        left_2 = parent_2[i][:cut_point]
        letters_left_2 = len(left_2.replace('-', ''))
        seen, cut_p1 = 0, 0
        
        if letters_left_2 > 0:
            for j, char in enumerate(parent_1[i]):
                if char != '-':
                    seen += 1
                if seen == letters_left_2:
                    cut_p1 = j + 1
                    break
        child_2.append(left_2 + parent_1[i][cut_p1:])
        
    m1 = max(len(s) for s in child_1)
    child_1 = [s.ljust(m1, '-') for s in child_1]
    m2 = max(len(s) for s in child_2)
    child_2 = [s.ljust(m2, '-') for s in child_2]
    
    return child_1, child_2


def apply_gap_mutation(alignment: list[str], mutation_rate: float = 0.05) -> list[str]:
    """Applies a gap shift mutation to introduce local structural variations."""
    mutated = copy.deepcopy(alignment)
    for i in range(len(mutated)):
        if random.random() < mutation_rate:
            seq = list(mutated[i])
            for _ in range(10):
                pos = random.randint(0, len(seq) - 2)
                if (seq[pos] == '-' and seq[pos+1] != '-') or (seq[pos] != '-' and seq[pos+1] == '-'):
                    seq[pos], seq[pos+1] = seq[pos+1], seq[pos]
                    break
            mutated[i] = "".join(seq)
    return mutated


def extract_qscore(terminal_output: str) -> tuple[str, str]:
    """Extracts Q-Score and TC-Score metrics from the QScore benchmark tool output."""
    match_q = re.search(r'Q=([0-9.]+)', terminal_output)
    match_tc = re.search(r'TC=([0-9.]+)', terminal_output)
    q = match_q.group(1) if match_q else "Error"
    tc = match_tc.group(1) if match_tc else "Error"
    return q, tc


def run_genetic_algorithm(fasta_path: str) -> tuple[str, float, str]:
    """
    Executes the genetic optimization pipeline for a single input instance:
    1. KAlign initialization (Warm Start).
    2. Population evolution over generations.
    3. Final alignment formatting and scoring.
    """
    base_dir = os.path.dirname(fasta_path)
    seq_name = os.path.basename(fasta_path).replace(".tfa", "").replace("_base.fasta", "")
    
    # Paths setup
    os.makedirs("base", exist_ok=True)
    os.makedirs("genetic-alg", exist_ok=True)
    
    temp_base_fasta = os.path.join("base", f"{seq_name}_base.fasta")
    output_ag_fasta = os.path.join("genetic-alg", f"{seq_name}_ag.fasta")
    
    # 1. KAlign Warm Start
    subprocess.run(["kalign", "-i", fasta_path, "-o", temp_base_fasta], capture_output=True)
    
    headers, base_seqs = read_aligned_fasta(temp_base_fasta)
    if not base_seqs:
        raise ValueError(f"Failed to read KAlign output for {seq_name}")
        
    num_seqs = len(base_seqs)
    length = len(base_seqs[0])
    
    generations = 100
    pop_size = 20
    
    log_buffer = []
    log_buffer.append(f"=== OPTIMIZATION LOG: {seq_name} ===\n")
    log_buffer.append(f"Dimensions: {num_seqs} sequences x {length} columns\n")
    log_buffer.append(f"Parameters: {generations} Generations | Population Size: {pop_size}\n\n")
    
    # Initialize Population with Alpha Individual (KAlign) + Mutated Clones
    population = [base_seqs] + [apply_gap_mutation(base_seqs, 0.40) for _ in range(pop_size - 1)]
    
    best_overall_score = float('-inf')
    best_overall_alignment = base_seqs
    
    for g in range(generations):
        pop_evaluated = [(calculate_wsp_fitness(ind), ind) for ind in population]
        pop_evaluated.sort(reverse=True, key=lambda x: x[0])
        
        current_best_score, current_best_ind = pop_evaluated[0]
        if current_best_score > best_overall_score:
            best_overall_score = current_best_score
            best_overall_alignment = current_best_ind
            
        new_population = [pop_evaluated[0][1], pop_evaluated[1][1]]
        
        if g % 10 == 0:
            log_msg = f"Generation {g:03d}/{generations} | Best WSP: {current_best_score:.2f}\n"
            log_buffer.append(log_msg)
            
        while len(new_population) < pop_size:
            half_pop = max(2, pop_size // 2)
            p1 = random.choice(pop_evaluated[:half_pop])[1]
            p2 = random.choice(pop_evaluated[:half_pop])[1]
            f1, f2 = conservative_crossover(p1, p2)
            new_population.extend([apply_gap_mutation(f1, 0.10), apply_gap_mutation(f2, 0.10)])
            
        population = new_population[:pop_size]
        
    log_buffer.append("\n=== GENETIC ALGORITHM COMPLETED ===\n")
    
    # Save optimized alignment
    save_fasta_result(headers, best_overall_alignment, output_ag_fasta)
    
    # Evaluate using reference benchmark tool (QScore)
    ref_path = f"bench/bench1.0/bali3/ref/{seq_name}"
    qscore_cmd = ["./qscore", "-seqdiffwarn", "-ignoretestcase", "-test", output_ag_fasta, "-ref", ref_path, "-truncname"]
    process = subprocess.run(qscore_cmd, capture_output=True, text=True)
    q_score, tc_score = extract_qscore(process.stdout)
    
    log_summary = "".join(log_buffer)
    formatted_output_alignment = "\n".join([f"{headers[i]}\n{best_overall_alignment[i]}" for i in range(len(best_overall_alignment))])
    
    return formatted_output_alignment, best_overall_score, log_summary