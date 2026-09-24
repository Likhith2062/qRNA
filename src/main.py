from dataclasses import dataclass
import time

from vienna import Vienna
from stem_generator import StemGenerator
from qubo_builder import QUBOBuilder
from quantum_solver import QuantumSolver
from decoder import Decoder
from metrics import evaluate


@dataclass
class PredictionResult:
    """All information produced by one RNA prediction run."""

    sequence: str
    stems: list
    qubo: dict
    prediction: str
    reference: str
    mfe: float
    prediction_energy: float
    metrics: dict
    prediction_matrix: object
    reference_matrix: object
    qubo_variables: int
    runtime: float


def predict(sequence: str) -> PredictionResult:
    """Run the complete RNA prediction pipeline."""
    start_time = time.perf_counter()
    sequence = sequence.strip().upper()

    if not sequence:
        raise ValueError("RNA sequence cannot be empty.")

    vienna = Vienna(sequence)

    generator = StemGenerator(vienna)
    stems = generator.generate()

    builder = QUBOBuilder(stems, vienna)
    Q = builder.build()

    solver = QuantumSolver(
        use_quantum_hardware=False,
        num_reads=100,
    )
    result = solver.solve(Q)

    decoder = Decoder(stems, len(sequence))
    prediction = decoder.dot_bracket(result.solution)
    prediction_matrix = decoder.adjacency_matrix(result.solution)

    prediction_energy = vienna.evaluate(prediction)

    reference, mfe = vienna.mfe()
    reference_matrix = Decoder.dot_bracket_to_matrix(reference)

    results = evaluate(reference_matrix, prediction_matrix)
    runtime = time.perf_counter() - start_time

    return PredictionResult(
        sequence=sequence,
        stems=stems,
        qubo=Q,
        prediction=prediction,
        reference=reference,
        mfe=mfe,
        prediction_energy=prediction_energy,
        metrics=results,
        prediction_matrix=prediction_matrix,
        reference_matrix=reference_matrix,
        qubo_variables=len(stems),
        runtime=runtime,
    )


def main():
    """Command-line interface."""
    sequence = input("RNA Sequence: ").strip().upper()

    try:
        result = predict(sequence)
    except Exception as exc:
        print(f"Error: {exc}")
        return

    print()
    print("Benchmark")
    print("---------")
    print(f"Sequence Length   : {len(result.sequence)}")
    print(f"Candidate Stems   : {len(result.stems)}")
    print(f"QUBO Variables    : {result.qubo_variables}")
    print(f"QUBO Coefficients : {len(result.qubo)}")
    print(f"Runtime           : {result.runtime:.4f} s")

    print()
    print("QUBO Predicted Structure")
    print("------------------------")
    print(result.prediction)
    print(f"Energy : {result.prediction_energy:.2f} kcal/mol")

    print()
    print("ViennaRNA MFE")
    print("-------------")
    print(result.reference)
    print(f"Energy : {result.mfe:.2f} kcal/mol")

    print()
    print("Evaluation")
    print("----------")

    for key, value in result.metrics.items():
        print(f"{key:>4} : {value}")


def launch_gui():
    """Launch the GUI without making PySide6 a CLI dependency."""
    try:
        from gui import run_gui
    except ImportError as exc:
        print(
            "GUI dependencies are not installed or gui.py is missing.\n"
            "Install PySide6 with: pip install PySide6"
        )
        print(f"Details: {exc}")
        return

    run_gui()


if __name__ == "__main__":
    main()
