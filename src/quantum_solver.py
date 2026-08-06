from dwave.system import DWaveSampler, EmbeddingComposite
import neal  # Simulated Annealing Sampler for local testing
from typing import Dict, List, Tuple

class QuantumRNAFoldSolver:
    def __init__(self, use_quantum_hardware: bool = False):
        self.use_hardware = use_quantum_hardware
        
        if self.use_hardware:
            # Requires valid DWAVE_API_TOKEN configured in environment
            self.sampler = EmbeddingComposite(DWaveSampler())
        else:
            # Local classical simulation of quantum annealing
            self.sampler = neal.SimulatedAnnealingSampler()

    def solve(self, Q: Dict[Tuple[int, int], float], num_reads: int = 100) -> List[int]:
        """
        Solves the QUBO problem and returns the best binary variable configuration.
        """
        response = self.sampler.sample_qubo(Q, num_reads=num_reads)
        best_sample = response.first.sample
        
        # Format the result as an ordered binary list
        max_idx = max(max(k) for k in Q.keys()) if Q else 0
        binary_solution = [best_sample.get(i, 0) for i in range(max_idx + 1)]
        
        return binary_solution
