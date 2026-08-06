import numpy as np
from typing import List, Tuple, Dict

class RNAQUBOBuilder:
    def __init__(self, sequence: str, stems: List[Tuple[int, int, int]]):
        """
        :param sequence: RNA sequence string (e.g., "ACGGUCAGU...")
        :param stems: List of stems represented as (i, j, length)
        """
        self.sequence = sequence
        self.stems = stems
        self.num_stems = len(stems)

    def _check_overlap(self, stem1: Tuple[int, int, int], stem2: Tuple[int, int, int]) -> bool:
        """Checks if two stems share the same nucleotide positions."""
        i1, j1, l1 = stem1
        i2, j2, l2 = stem2
        
        indices1 = set(range(i1, i1 + l1)).union(set(range(j1 - l1 + 1, j1 + 1)))
        indices2 = set(range(i2, i2 + l2)).union(set(range(j2 - l2 + 1, j2 + 1)))
        
        return len(indices1.intersection(indices2)) > 0

    def _check_pseudoknot(self, stem1: Tuple[int, int, int], stem2: Tuple[int, int, int]) -> bool:
        """Checks if two stems cross each other forming a pseudoknot."""
        i1, j1, _ = stem1
        i2, j2, _ = stem2
        
        return (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1)

    def build_m1_qubo(self, alpha: float = 1.0, beta: float = 0.5, P_overlap: float = 10.0, P_pk: float = 2.0) -> Dict[Tuple[int, int], float]:
        """
        Model 1 (M1 - Baseline Stem-level):
        Maximizes stem length and base pairs while penalizing overlap and pseudoknots.
        """
        Q = {}

        # Linear terms (stem energy/length benefits)
        for idx, (i, j, length) in enumerate(self.stems):
            # Reward longer stems
            reward = -(alpha * length + beta)
            Q[(idx, idx)] = reward

        # Quadratic terms (conflict penalties)
        for idx1 in range(self.num_stems):
            for idx2 in range(idx1 + 1, self.num_stems):
                stem1 = self.stems[idx1]
                stem2 = self.stems[idx2]

                penalty = 0.0
                if self._check_overlap(stem1, stem2):
                    penalty += P_overlap
                elif self._check_pseudoknot(stem1, stem2):
                    penalty += P_pk

                if penalty > 0:
                    Q[(idx1, idx2)] = penalty

        return Q

    def build_m3_qubo(self, weights: Dict[str, float]) -> Dict[Tuple[int, int], float]:
        """
        Model 3 (M3 - Physics-inspired QUBO):
        Uses thermodynamics weights, hairpin penalties, and polymer physics-based pseudoknot penalties.
        """
        w_energy = weights.get("w_energy", 1.0)
        P_overlap = weights.get("P_overlap", 15.0)
        P_pk = weights.get("P_pk", 1.5)

        Q = {}

        # Linear energy terms (pseudo Nearest-Neighbor energy approximation)
        for idx, (i, j, length) in enumerate(self.stems):
            # Simplification: longer stems yield lower (better) free energy
            approx_free_energy = -1.2 * length
            Q[(idx, idx)] = w_energy * approx_free_energy

        # Quadratic terms
        for idx1 in range(self.num_stems):
            for idx2 in range(idx1 + 1, self.num_stems):
                stem1 = self.stems[idx1]
                stem2 = self.stems[idx2]

                if self._check_overlap(stem1, stem2):
                    Q[(idx1, idx2)] = P_overlap
                elif self._check_pseudoknot(stem1, stem2):
                    # Distance-dependent pseudoknot penalty approximation
                    dist = abs(stem1[0] - stem2[0])
                    Q[(idx1, idx2)] = P_pk * (1.0 + 0.1 * np.log(dist + 1))

        return Q
