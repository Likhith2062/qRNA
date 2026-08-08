from vienna import Vienna
from stem_generator import StemGenerator
from qubo_builder import QUBOBuilder
from quantum_solver import QuantumSolver
from decoder import Decoder
from metrics import evaluate


def main():

    # -----------------------------------------------------
    # Input sequence
    # -----------------------------------------------------

    sequence = input(
        "RNA Sequence: "
    ).strip().upper()

    vienna = Vienna(sequence)

    # -----------------------------------------------------
    # Candidate stems
    # -----------------------------------------------------

    generator = StemGenerator(vienna)

    stems = generator.generate()

    print(f"\nGenerated {len(stems)} candidate stems.")

    # -----------------------------------------------------
    # Build QUBO
    # -----------------------------------------------------

    builder = QUBOBuilder(
        stems,
        vienna,
    )

    Q = builder.build()

    print(
        f"QUBO contains {len(Q)} coefficients."
    )

    # -----------------------------------------------------
    # Solve
    # -----------------------------------------------------

    solver = QuantumSolver(
        use_quantum_hardware=False,
        num_reads=100,
    )

    result = solver.solve(Q)

    print(
        f"Best QUBO energy : {result.energy:.4f}"
    )

    # -----------------------------------------------------
    # Decode
    # -----------------------------------------------------

    decoder = Decoder(
        stems,
        len(sequence),
    )

    prediction = decoder.dot_bracket(
        result.solution
    )

    prediction_matrix = decoder.adjacency_matrix(
        result.solution
    )

    print()

    print("Predicted Structure")

    print(prediction)

    # -----------------------------------------------------
    # ViennaRNA reference
    # -----------------------------------------------------

    reference, mfe = vienna.mfe()

    reference_matrix = Decoder.dot_bracket_to_matrix(
        reference
    )

    print()

    print("ViennaRNA MFE")

    print(reference)

    print(
        f"Energy : {mfe:.2f} kcal/mol"
    )

    # -----------------------------------------------------
    # Evaluation
    # -----------------------------------------------------

    results = evaluate(
        reference_matrix,
        prediction_matrix,
    )

    print()

    print("Evaluation")

    print("--------------------")

    for key, value in results.items():

        print(
            f"{key:>4} : {value}"
        )


if __name__ == "__main__":

    main()