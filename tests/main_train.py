import numpy as np
from src.preprocessing import find_possible_stems
from src.qubo_models import RNAQUBOBuilder
from src.quantum_solver import QuantumRNAFoldSolver
from src.metrics import calculate_mcc

def main():
    # Sample RNA Sequence
    rna_sequence = "ACGGUCAGUCCUUUACUGA"
    print(f"[1/4] Processing Sequence: {rna_sequence}")

    # 1. Extract Candidate Stems
    stems = find_possible_stems(rna_sequence, min_stem_length=2, min_loop_length=3)
    print(f"[2/4] Found {len(stems)} candidate stems.")

    # 2. Build QUBO Matrix
    qubo_builder = RNAQUBOBuilder(rna_sequence, stems)
    weights = {"w_energy": 1.0, "P_overlap": 20.0, "P_pk": 2.0}
    Q_matrix = qubo_builder.build_m3_qubo(weights)

    # 3. Solve QUBO (using local simulated annealer)
    print("[3/4] Solving QUBO via Simulated Quantum Annealing...")
    solver = QuantumRNAFoldSolver(use_quantum_hardware=False)
    solution = solver.solve(Q_matrix, num_reads=50)

    # 4. Filter selected stems
    selected_stems = [stems[i] for i, selected in enumerate(solution) if selected == 1]
    print(f"[4/4] Optimization Complete. Selected {len(selected_stems)} stems:")
    for stem in selected_stems:
        print(f"  - Stem starting at i={stem[0]}, j={stem[1]} (length {stem[2]})")

if __name__ == "__main__":
    main()
