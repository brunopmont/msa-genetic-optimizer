# =====================================================================
# MSA Genetic Optimizer - Makefile
# =====================================================================

PYTHON = python3

.PHONY: all base benchmark clean help

all: benchmark

help:
	@echo "Available commands:"
	@echo "  make base       - Generate initial seed alignments using KAlign"
	@echo "  make benchmark  - Run the complete Genetic Algorithm benchmark suite"
	@echo "  make clean      - Remove temporary files, logs, and build artifacts"

base:
	@echo "Generating base alignments..."
	$(PYTHON) gerar_base.py

benchmark:
	@echo "Running full benchmark suite..."
	$(PYTHON) benchmark_completo.py

clean:
	@echo "Cleaning up temporary directories..."
	rm -rf __pycache__
	rm -rf base/
	rm -rf genetic-alg/
	rm -rf logs_individuais/
	rm -rf resultados/
	@echo "Cleanup completed successfully."