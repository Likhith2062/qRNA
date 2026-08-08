from typing import Dict, List, Tuple
from dataclasses import dataclass

@dataclass(slots=True)
class SolverResult:
    """
    Result returned by a QUBO solver.
    """

    solution: list[int]

    energy: float

class QuantumSolver:
    def __init__(
        self,
        use_quantum_hardware: bool = False,
        num_reads: int = 100,
    ):
        """
        Parameters
        ----------
        use_quantum_hardware
            If True, use a D-Wave quantum annealer.

        num_reads
            Number of annealing runs.
        """

        self.use_hardware = use_quantum_hardware

        self.num_reads = num_reads

        if self.use_hardware:
            from dwave.system import (
                DWaveSampler,
                EmbeddingComposite,
            )

            self.sampler = EmbeddingComposite(
                DWaveSampler()
            )

        else:
            try:
                from dwave.samplers import SimulatedAnnealingSampler
                self.sampler = SimulatedAnnealingSampler()
            except ImportError as exc:
                raise ImportError(
                    "Local simulated annealer not found.\n"
                    "Install it with:\n"
                    "pip install dwave-neal\n"
                    "or\n"
                    "pip install dwave-ocean-sdk"
                ) from exc

    def solve(self, Q: Dict[Tuple[int, int], float]) -> SolverResult:
        """
        Solves the QUBO problem and returns the best binary variable configuration.
        """
        response = self.sampler.sample_qubo(Q, num_reads=self.num_reads)
        best_sample = response.first.sample
        
        # Format the result as an ordered binary list
        max_idx = max(max(k) for k in Q.keys()) if Q else 0
        binary_solution = [best_sample.get(i, 0) for i in range(max_idx + 1)]
        
        return SolverResult(solution=binary_solution, energy=response.first.energy)
