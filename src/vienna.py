"""
vienna.py
=========

Thin wrapper around ViennaRNA.

Responsibilities
----------------
- Hold a FoldCompound object
- Expose thermodynamic lookup tables
- Evaluate RNA secondary structures
- Evaluate local motifs
- Draw RNA structures using ViennaRNA's native plotting routines
"""

from __future__ import annotations

from pathlib import Path

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
        """Convert ViennaRNA energy units to kcal/mol."""
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

    def mfe(self) -> tuple[str, float]:
        """
        Compute the minimum free-energy (MFE) secondary structure.

        Returns
        -------
        tuple[str, float]
            (dot_bracket, energy)
        """

        structure, energy = self.fc.mfe()

        return structure, float(energy)

    # ---------------------------------------------------------
    # Local motif energies
    # ---------------------------------------------------------

    def stack_energy(self, i: int, j: int) -> float:

        return self._to_kcal(self.fc.E_stack(i, j))

    def stem_stacking_energy(self, stem) -> float:
        """
        Return the total nearest-neighbor stacking energy
        of a candidate stem.
        """

        energy = 0.0

        for offset in range(stem.length - 1):

            i = stem.start5 + offset
            j = stem.end3 - offset

            e = self.stack_energy(i, j)

            if e < 1e6:
                energy += e

        return energy

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

    def plot(self, dot_bracket: str) -> str:
        """
        Generate a ViennaRNA structure drawing entirely in memory.

        Parameters
        ----------
        dot_bracket : str
            RNA secondary structure in dot-bracket notation.

        Returns
        -------
        str
            SVG document containing the ViennaRNA-rendered structure.

        Notes
        -----
        Nothing is written to disk by this method. The returned SVG string
        can be displayed directly by the GUI or saved later when the user
        explicitly chooses a location.
        """

        if len(dot_bracket) != len(self.sequence):
            raise ValueError(
                "Dot-bracket structure length must match the RNA sequence."
            )

        # ViennaRNA's Python binding writes the SVG to a file and
        # returns a status code rather than the SVG contents. Use a
        # temporary file only as an internal bridge, immediately read
        # the SVG into memory, and delete the file.
        import os
        import tempfile

        fd, filename = tempfile.mkstemp(suffix=".svg")
        os.close(fd)

        try:
            result = RNA.svg_rna_plot(
                self.sequence,
                dot_bracket,
                filename,
            )

            if not os.path.exists(filename):
                raise RuntimeError(
                    f"ViennaRNA failed to create the SVG file "
                    f"(status: {result})."
                )

            with open(filename, "r", encoding="utf-8") as file:
                svg = file.read()

            if not svg.strip():
                raise RuntimeError(
                    "ViennaRNA generated an empty SVG document."
                )

            return svg

        finally:
            # The SVG is now safely held in memory. Remove the temporary
            # file so no structure image is left on disk.
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

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
