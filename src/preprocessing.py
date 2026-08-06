import numpy as np
from typing import List, Tuple

# Paires de bases canoniques (Watson-Crick + Wobble G-U)
VALID_PAIRS = {('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C'), ('G', 'U'), ('U', 'G')}

def is_valid_pair(base1: str, base2: str) -> bool:
    """Vérifie si deux bases peuvent former une paire valide."""
    return (base1.upper(), base2.upper()) in VALID_PAIRS

def find_possible_stems(sequence: str, min_stem_length: int = 2, min_loop_length: int = 3) -> List[Tuple[int, int, int]]:
    """
    Identifie toutes les tiges (stems) potentielles dans une séquence d'ARN.
    
    Retourne une liste de tuples : (index_i, index_j, longueur_tige)
    où i est le début du brin 5' et j est le début du brin 3'.
    """
    n = len(sequence)
    stems = []

    for i in range(n):
        for j in range(i + min_loop_length + 2 * min_stem_length - 1, n):
            # Vérifier la longueur de la tige contiguë
            length = 0
            while (i + length < j - length - min_loop_length) and is_valid_pair(sequence[i + length], sequence[j - length]):
                length += 1

            if length >= min_stem_length:
                stems.append((i, j, length))

    return stems

if __name__ == "__main__":
    # Séquence ARN de test
    sample_rna = "ACGGUCAGUCCUUUACUGA"
    print(f"Analyse de la séquence ARN : {sample_rna}")
    
    stems = find_possible_stems(sample_rna)
    print(f"Tiges détectées (i, j, longueur) : {stems}")
