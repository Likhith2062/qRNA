from dataclasses import dataclass

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
    solver_energy: float
    metrics: dict
    prediction_matrix: object
    reference_matrix: object


def predict(sequence: str) -> PredictionResult:
    """Run the complete RNA prediction pipeline."""
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

    reference, mfe = vienna.mfe()
    reference_matrix = Decoder.dot_bracket_to_matrix(reference)

    results = evaluate(reference_matrix, prediction_matrix)

    return PredictionResult(
        sequence=sequence,
        stems=stems,
        qubo=Q,
        prediction=prediction,
        reference=reference,
        mfe=mfe,
        solver_energy=result.energy,
        metrics=results,
        prediction_matrix=prediction_matrix,
        reference_matrix=reference_matrix,
    )


def main():
    """Command-line interface."""
    sequence = input("RNA Sequence: ").strip().upper()

    try:
        result = predict(sequence)
    except Exception as exc:
        print(f"Error: {exc}")
        return

    print(f"\nGenerated {len(result.stems)} candidate stems.")
    print(f"QUBO contains {len(result.qubo)} coefficients.")

    print()
    print("QUBO Predicted Structure")
    print("------------------------")
    print(result.prediction)

    print()
    print("ViennaRNA MFE")
    print("-------------")
    print(result.reference)
    print(f"Energy : {result.mfe:.2f} kcal/mol")

    print()
    print("Solver")
    print("------")
    print(f"QUBO Energy : {result.solver_energy:.4f}")

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

    run_gui(predict)


if __name__ == "__main__":
    main()
