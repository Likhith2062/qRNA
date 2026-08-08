"""
decoder.py
==========

Decode a QUBO solution into RNA secondary structure
representations.

The decoder performs no optimisation.

It simply interprets the binary solution returned by
the quantum solver.
"""

from __future__ import annotations
import numpy as np
from stem_generator import Stem


class Decoder:

    def __init__(
        self,
        stems: list[Stem],
        sequence_length: int,
    ):

        self.stems = stems

        self.n = sequence_length

    def selected_stems(
        self,
        solution: list[int],
    ) -> list[Stem]:
        """
        Return every selected stem.

        solution[i] == 1

        means stem i is selected.
        """

        selected = []

        for bit, stem in zip(
            solution,
            self.stems,
        ):
            if bit:
                selected.append(
                    stem
                )

        return selected

    def base_pairs(
        self,
        solution: list[int],
    ) -> list[tuple[int, int]]:
        """
        Return every base pair represented by the
        selected stems.

        Returned coordinates are 1-based.
        """

        pairs = []

        for stem in self.selected_stems(
            solution
        ):

            for offset in range(
                stem.length
            ):

                left = (
                    stem.start5
                    +
                    offset
                )

                right = (
                    stem.end3
                    -
                    offset
                )

                pairs.append(
                    (
                        left,
                        right,
                    )
                )

        return sorted(pairs)

    def dot_bracket(
        self,
        solution: list[int],
    ) -> str:
        """
        Convert selected base pairs into
        dot-bracket notation.

        Currently supports only
        non-pseudoknotted structures.
        """

        structure = [

            "."

            for _ in range(self.n)

        ]

        for left, right in self.base_pairs(
            solution
        ):

            structure[left - 1] = "("

            structure[right - 1] = ")"

        return "".join(structure)

    def adjacency_matrix(
        self,
        solution: list[int],
    ):
        """
        Return the base-pair adjacency matrix.
        """

        matrix = [

            [0] * self.n

            for _ in range(self.n)

        ]

        for left, right in self.base_pairs(
            solution
        ):

            matrix[left - 1][right - 1] = 1

            matrix[right - 1][left - 1] = 1

        return np.array(
            matrix,
            dtype=np.int8,
        )