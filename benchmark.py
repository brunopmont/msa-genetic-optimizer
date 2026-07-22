"""
complete_benchmark.py
Automated benchmark runner for the Genetic Algorithm MSA optimizer 
across BAliBase 3 instances.
"""

import time
import os
from refinement_genetic_algorithm import run_genetic_algorithm

# Definition of the evaluation instances and their configurations
BENCHMARK_INSTANCES = [
    {"id": "BB11001", "path": "qscore/BB11001.tfa"},
    {"id": "BB11010", "path": "qscore/BB11010.tfa"},
    {"id": "BB12012", "path": "qscore/BB12012.tfa"},
    {"id": "BB12022", "path": "qscore/BB12022.tfa"},
    {"id": "BB20030", "path": "qscore/BB20030.tfa"},
    {"id": "BB30016", "path": "qscore/BB30016.tfa"},
    {"id": "BB20010", "path": "qscore/BB20010.tfa"},
    {"id": "BB30024", "path": "qscore/BB30024.tfa"},
    {"id": "BB40011", "path": "qscore/BB40011.tfa"},
    {"id": "BB50014", "path": "qscore/BB50014.tfa"},
]


def save_execution_log(instance_id: str, log_data: str) -> None:
    """Saves individual execution logs to the logs directory."""
    os.makedirs("logs_individuais", exist_ok=True)
    file_path = os.path.join("logs_individuais", f"log_execucao_{instance_id}.txt")
    with open(file_path, "w", encoding="utf-8") as file:
        file.write(log_data)


def run_benchmark_suite() -> None:
    """Executes the genetic algorithm across all benchmark instances sequentially."""
    os.makedirs("resultados", exist_ok=True)
    summary_path = os.path.join("resultados", "resultados_gerais_benchmark.txt")
    
    with open(summary_path, "w", encoding="utf-8") as summary_file:
        summary_file.write("Instance ID\tExecution Time (s)\tBest WSP Score\n")
        
        for instance in BENCHMARK_INSTANCES:
            instance_id = instance["id"]
            file_path = instance["path"]
            
            print(f"\n[INFO] Starting execution for instance: {instance_id}")
            start_time = time.time()
            
            try:
                # Runs the genetic algorithm optimization
                best_alignment, best_score, execution_log = run_genetic_algorithm(file_path)
                elapsed_time = time.time() - start_time
                
                print(f"[SUCCESS] {instance_id} completed in {elapsed_time:.2f}s | WSP: {best_score}")
                
                # Save individual results and logs
                save_execution_log(instance_id, execution_log)
                
                result_filename = os.path.join("resultados", f"resultado_final_{instance_id}.txt")
                with open(result_filename, "w", encoding="utf-8") as res_file:
                    res_file.write(best_alignment)
                
                summary_file.write(f"{instance_id}\t{elapsed_time:.2f}\t{best_score}\n")
                
            except Exception as e:
                print(f"[ERROR] Failed to process instance {instance_id}: {e}")
                summary_file.write(f"{instance_id}\tERROR\tERROR\n")


if __name__ == "__main__":
    print("=== MSA Genetic Optimizer - Benchmark Suite ===")
    run_benchmark_suite()
    print("\n[INFO] Benchmark suite execution finished.")