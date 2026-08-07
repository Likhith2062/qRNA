"""
vienna.py
==========

Thin wrapper around ViennaRNA.

Responsibilities
----------------
- Hold a FoldCompound object
- Expose thermodynamic lookup tables
- Evaluate RNA secondary structures
- Evaluate local motifs
- Draw RNA structures

This module intentionally DOES NOT expose any ViennaRNA
optimization routines (mfe, subopt, backtrack, etc.).
"""

from __future__ import annotations

import RNA


class Vienna:

    def __init__(self, sequence: str):

        if not sequence:
            raise ValueError("Sequence cannot be empty.")

        sequence = sequence.upper().replace("T", "U")

        allowed = {"A", "C", "G", "U"}

        invalid = set(sequence) - allowed

        if invalid:
            raise ValueError(
                f"Invalid nucleotide(s): {sorted(invalid)}"
            )

        self.sequence = sequence

        # ViennaRNA objects
        self.fc = RNA.fold_compound(sequence)
        self.params = self.fc.params
        self.model = self.params.model_details

        # Thermodynamic tables
        self.stack = self.params.stack
        self.hairpin = self.params.hairpin
        self.bulge = self.params.bulge
        self.internal_loop = self.params.internal_loop
        self.dangle5 = self.params.dangle5
        self.dangle3 = self.params.dangle3
        self.ninio = self.params.ninio

        # Pair lookup
        self.pair_table = self.model.pair
        self.reverse_pair = self.model.rtype

        # Model information
        self.temperature = self.params.temperature
        self.min_loop_size = self.model.min_loop_size

    # ---------------------------------------------------------
    # Utility
    # ---------------------------------------------------------

    @staticmethod
    def _to_kcal(value: int | float) -> float:
        """
        Convert ViennaRNA energy units to kcal/mol.
        """
        return float(value) / 100.0

    @staticmethod
    def encode_base(base: str) -> int:

        lookup = {
            "A": 1,
            "C": 2,
            "G": 3,
            "U": 4,
        }

        return lookup.get(base.upper(), 0)

    # ---------------------------------------------------------
    # Pair utilities
    # ---------------------------------------------------------

    def pair_type(self, base1: str, base2: str) -> int:

        a = self.encode_base(base1)
        b = self.encode_base(base2)

        return self.pair_table[a][b]

    def can_pair(self, base1: str, base2: str) -> bool:

        return self.pair_type(base1, base2) != 0

    # ---------------------------------------------------------
    # Whole structure
    # ---------------------------------------------------------

    def evaluate(self, dot_bracket: str) -> float:

        return float(self.fc.eval_structure(dot_bracket))

    # ---------------------------------------------------------
    # Local motif energies
    # ---------------------------------------------------------

    def stack_energy(self, i: int, j: int) -> float:

        return self._to_kcal(self.fc.E_stack(i, j))

    def hairpin_energy(self, i: int, j: int) -> float:

        return self._to_kcal(
            self.fc.eval_hp_loop(i, j)
        )

    def internal_loop_energy(
        self,
        i: int,
        j: int,
        k: int,
        l: int,
    ) -> float:

        return self._to_kcal(
            self.fc.eval_int_loop(i, j, k, l)
        )

    # ---------------------------------------------------------
    # Plot
    # ---------------------------------------------------------

    def plot(self, dot_bracket: str, filename="rna.ps"):

        RNA.file_PS_rnaplot(
            self.sequence,
            dot_bracket,
            filename,
        )

    # ---------------------------------------------------------
    # Information
    # ---------------------------------------------------------

    def info(self):

        return {
            "sequence": self.sequence,
            "length": len(self.sequence),
            "temperature": self.temperature,
            "min_loop_size": self.min_loop_size,
        }