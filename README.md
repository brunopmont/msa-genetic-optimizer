# MSA Genetic Optimizer

Operations Research final project exploring Multiple Sequence Alignment (MSA) optimization using a Genetic Algorithm with a deterministic warm start (KAlign) and conservative reproduction operators. Developed for the Computer Science Institute at Universidade Federal Fluminense (UFF).

## Overview

MSA is an NP-hard problem characterized by a combinatorial explosion of gap insertions across biological sequences. This project implements a hybrid Genetic Algorithm designed around Holland's Building Block Theory, prioritizing structural preservation over destructive stochastic perturbations commonly found in literature (such as Ant Colony Optimization with chaotic jumps).

- **Core Heuristic:** Warm start initialized via KAlign.
- **Objective Function:** Weighted Sum of Pairs (WSP) with affine gap penalties.
- **Benchmark:** Validated against BAliBase 3 across 10 distinct complexity tiers.

## Repository Structure

- `ag_refinamento.py`: Main implementation of the Genetic Algorithm, selection, crossover, and mutation operators.
- `benchmark_completo.py`: Automated benchmarking script for running instances sequentially.
- `gerar_base.py` / `gerar_amostra.py`: Data preparation and formatting utilities.
- `data_analysis.ipynb`: Analysis notebook and performance plotting scripts.
- `resultados/`: Output text files and convergence plots.
- `logs_individuais/`: Raw execution logs per BAliBase instance.

## Requirements

- Python 3.10+
- Dependencies: Standard scientific libraries (Pandas, NumPy, Matplotlib, Seaborn).
- External Tool: `kalign` and the `qscore` benchmark evaluation binary.

## Usage

To generate the initial seed alignments and execute the complete benchmark suite locally:

```bash
make base
make benchmark
make clean
```

## Results summary

The proposed GA matches or outperforms the reference ACO literature in 80% of the tested BAliBase instances while consuming only 1/3 of the computational budget (2,000 FEs vs. 6,000 FEs).