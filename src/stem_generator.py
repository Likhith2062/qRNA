"""
stem_generator.py
=================

Generate every possible candidate stem from an RNA sequence.

A stem is purely geometric.
No thermodynamic calculations are performed here.

Each generated Stem corresponds to ONE binary variable
in the later QUBO formulation.
"""

from __future__ import annotations

from dataclasses import dataclass
from vienna import Vienna


# ---------------------------------------------------------
# Data structure
# ---------------------------------------------------------

@dataclass(slots=True)
class Stem:
    """
    Candidate stem.

    Coordinates are 1-based to match ViennaRNA.
    """

    id: int

    # 5' strand
    start5: int
    end5: int

    # 3' strand
    start3: int
    end3: int

    # consecutive base pairs
    length: int


# ---------------------------------------------------------
# Stem Generator
# ---------------------------------------------------------

class StemGenerator:

    def __init__(
        self,
        vienna: Vienna,
        min_stem_length: int = 2,
    ):

        self.vienna = vienna

        self.sequence = vienna.sequence

        self.n = len(self.sequence)

        self.min_stem = min_stem_length

        self.stems = []

        self.next_id = 0

    # -----------------------------------------------------

    def generate(self):

        """
        Generate every valid stem.
        """

        self.stems.clear()

        self.next_id = 0

        for i in range(self.n):

            for j in range(i + self.vienna.min_loop_size + 1,
                           self.n):

                if not self.vienna.can_pair(
                    self.sequence[i],
                    self.sequence[j]
                ):
                    continue

                self._scan_diagonal(i, j)

        return self.stems

    # -----------------------------------------------------

    def _scan_diagonal(self, i, j):

        """
        Extend a possible stem diagonally.
        """

        pairs = []

        left = i
        right = j

        while (
            left < right
            and (right - left - 1) >= self.vienna.min_loop_size
            and self.vienna.can_pair(
                self.sequence[left],
                self.sequence[right]
            )
        ):

            pairs.append((left + 1, right + 1))

            left += 1
            right -= 1

        if len(pairs) < self.min_stem:
            return

        self._generate_substems(pairs)

    # -----------------------------------------------------

    def _generate_substems(self, pairs):

        """
        Generate every contiguous substem.

        Example

        [(1,10),(2,9),(3,8),(4,7)]

        gives

        length4

        length3
        length3

        length2
        length2
        length2
        """

        L = len(pairs)

        for length in range(
            self.min_stem,
            L + 1,
        ):

            for start in range(
                L - length + 1
            ):

                sub = pairs[
                    start:
                    start + length
                ]

                first = sub[0]

                last = sub[-1]

                stem = Stem(

                    id=self.next_id,

                    start5=first[0],
                    end5=last[0],

                    start3=last[1],
                    end3=first[1],

                    length=length,
                )

                self.stems.append(stem)

                self.next_id += 1