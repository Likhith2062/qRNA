"""
qubo_builder.py
===============

Construct the QUBO matrix for RNA secondary structure prediction.

This implementation follows Model 3 of

Zaborniak et al.
Quantum Annealing for RNA Secondary Structure Prediction.

Current implementation
----------------------

✓ Linear Hamiltonian

Future commits

- Overlap penalty
- Pseudoknot penalty
"""

from __future__ import annotations

from dataclasses import dataclass

from vienna import Vienna
from stem_generator import Stem

import math

# ---------------------------------------------------------
# Parameters
# ---------------------------------------------------------

@dataclass(frozen=True)
class QUBOParameters:
    """
    Model 3 coefficients.

    Values reported after SPSA optimisation
    in the published paper.
    """

    alpha: float = 1.604

    beta: float = 2.212

    p1: float = 1.495

    p2: float = 1.338

    overlap_penalty: float = 1000.0

    nucleotide_length: float = 6.5

# ---------------------------------------------------------
# Cached stem data
# ---------------------------------------------------------

@dataclass(slots=True)
class StemData:
    """
    Cached information for one candidate stem.

    Thermodynamic quantities are computed exactly once
    during construction of the QUBO builder.
    """

    stem: Stem

    stacking_energy: float

    hairpin_penalty: float

@dataclass(slots=True)
class PseudoknotGeometry:
    """
    Geometric description of an H-type pseudoknot.
    """

    stem1: Stem

    stem2: Stem

    gap5: int

    gap3: int

    single_stranded: int

# ---------------------------------------------------------
# Builder
# ---------------------------------------------------------

class QUBOBuilder:

    # -------------------------------------------------
    # Hajdin Table S1
    # -------------------------------------------------

    LAMBDA_IL = {

        2: 0,

        3: 24,

        4: 24,

        5: 34,

        6: 24,

        7: 35,

        8: 58,

        9: 82,

        10: 65,

        11: 527,

        12: 2447,

        13: 4199,

        14: 6564,

        15: 12540,
    }

    # Penalty for pseudoknots outside the
    # experimentally calibrated region.

    INVALID_PSEUDOKNOT = 1e12

    def __init__(
        self,
        stems: list[Stem],
        vienna: Vienna,
        parameters: QUBOParameters = QUBOParameters(),
    ):

        self.stems = stems

        self.vienna = vienna

        self.params = parameters

        self.Q: dict[tuple[int, int], float] = {}

        # -------------------------------------------------
        # Cache thermodynamic data
        # -------------------------------------------------

        self.data: dict[int, StemData] = {}

        for stem in self.stems:

            self.data[stem.id] = StemData(

                stem=stem,

                stacking_energy=
                    self.vienna.stem_stacking_energy(stem),

                hairpin_penalty=
                    self.vienna.hairpin_energy(
                        stem.end5,
                        stem.start3,
                    ),
            )

        # Largest stem energy (paper notation μ)

        self.mu = self._compute_mu()

    # -------------------------------------------------

    def _compute_mu(self) -> float:
        """
        Maximum stem stacking energy.
        """

        if not self.data:
            return 0.0

        return max(

            -data.stacking_energy

            for data in self.data.values()

        )

    # -------------------------------------------------

    def _linear_term(self, stem: Stem) -> float:
        """
        Equation (8)

        α(k−μ)^2 − β(k−l)
        """
        data = self.data[stem.id]

        # Convert thermodynamic free energy (negative)
        # into stabilizing energy magnitude (positive),
        # as assumed by the Model 3 Hamiltonian.

        # ViennaRNA reports stabilizing stem free energies
        # as negative ΔG values.
        #
        # Model 3 uses stem energy as a positive stabilizing
        # quantity. Therefore convert
        #
        #     k = -ΔG
        #
        # before evaluating the Hamiltonian.

        k = -data.stacking_energy

        l = data.hairpin_penalty

        return (
            self.params.alpha * (k - self.mu) ** 2
            - self.params.beta * (k - l)
        )

    # -------------------------------------------------

    def build(self) -> dict[tuple[int, int], float]:
        """
        Construct the linear QUBO.

        Quadratic terms are added
        in later commits.
        """

        self.Q.clear()

        for stem in self.stems:

            self.Q[(stem.id, stem.id)] = (
                self._linear_term(stem)
            )

        self._quadratic_overlap()

        self._quadratic_pseudoknots()

        return self.Q

    def _stems_overlap(
        self,
        stem1: Stem,
        stem2: Stem,
    ) -> bool:
        """
        Return True if two stems share one or more nucleotides.

        Overlapping stems cannot coexist in the final structure.
        """

        used1 = set()

        for i in range(stem1.length):

            used1.add(stem1.start5 + i)

            used1.add(stem1.end3 - i)

        for i in range(stem2.length):

            if (stem2.start5 + i) in used1:
                return True

            if (stem2.end3 - i) in used1:
                return True

        return False

    def _pseudoknot_geometry(
        self,
        stem1: Stem,
        stem2: Stem,
    ) -> PseudoknotGeometry | None:
        """
        Return the geometry of an H-type pseudoknot.

        If the stems do not form a pseudoknot,
        return None.
        """

        #
        # Normalize ordering
        #

        if stem2.start5 < stem1.start5:
            stem1, stem2 = stem2, stem1

        #
        # Crossing condition
        #

        if not (
            stem1.start5
            < stem2.start5
            < stem1.start3
            < stem2.start3
        ):
            return None

        #
        # Single-stranded gaps
        #

        gap5 = stem2.start5 - stem1.end5 - 1

        gap3 = stem2.end3 - stem1.start3 - 1

        #
        # Invalid geometry
        #

        if gap5 < 0 or gap3 < 0:
            return None

        return PseudoknotGeometry(

            stem1=stem1,

            stem2=stem2,

            gap5=gap5,

            gap3=gap3,

            single_stranded=gap5 + gap3,
        )
    
    def _lambda(
        self,
        stem: Stem,
    ) -> float:
        """
        Return the Hajdin Table S1 penalty constant
        for one in-line helix.
        """

        if stem.length > 15:
            return self.INVALID_PSEUDOKNOT

        return self.LAMBDA_IL[stem.length]

    def _pseudoknot_penalty(
        self,
        geometry: PseudoknotGeometry,
    ) -> float:
        """
        Polymer entropy penalty from Equation (7).
        """

        lambda1 = self._lambda(
            geometry.stem1
        )

        lambda2 = self._lambda(
            geometry.stem2
        )

        #
        # Unsupported pseudoknot
        #

        if (
            lambda1 == self.INVALID_PSEUDOKNOT
            or
            lambda2 == self.INVALID_PSEUDOKNOT
        ):

            return self.INVALID_PSEUDOKNOT

        #
        # Prevent log(0)
        #

        if geometry.single_stranded <= 0:
            return self.INVALID_PSEUDOKNOT

        #
        # Prevent log(0) for lambda term
        #

        lambda_sum = lambda1 + lambda2

        if lambda_sum <= 0:

            return self.INVALID_PSEUDOKNOT

        entropy = (

            self.params.p1

            *

            (
                2.0
                +
                math.log(
                    geometry.single_stranded
                )
            )

            +

            self.params.p2

            *

            math.log(
                lambda_sum
            )

        )

        return entropy

    def _quadratic_overlap(self):

        for i in range(len(self.stems)):

            stem1 = self.stems[i]

            for j in range(i + 1, len(self.stems)):

                stem2 = self.stems[j]

                if self._stems_overlap(stem1, stem2):

                    self.Q[(stem1.id, stem2.id)] = (
                        self.params.overlap_penalty
                    )

    def _quadratic_pseudoknots(self):
        """
        Add entropy-based pseudoknot penalties.
        """

        for i in range(len(self.stems)):

            stem1 = self.stems[i]

            for j in range(i + 1, len(self.stems)):

                stem2 = self.stems[j]

                geometry = self._pseudoknot_geometry(
                    stem1,
                    stem2,
                )

                if geometry is None:
                    continue

                penalty = self._pseudoknot_penalty(
                    geometry
                )

                #
                # Ignore unsupported pseudoknots
                #

                if penalty == self.INVALID_PSEUDOKNOT:
                    continue

                #
                # If overlap already exists,
                # accumulate the penalties.
                #

                key = (
                    geometry.stem1.id,
                    geometry.stem2.id,
                )

                self.Q[key] = (

                    self.Q.get(key, 0.0)

                    +

                    penalty

                )