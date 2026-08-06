# RNA Folding using Variational Hybrid Quantum Annealing (QUBO)

This project implements and compares various QUBO (*Quadratic Unconstrained Binary Optimization*) formulation models applied to the RNA secondary structure folding problem including pseudoknots. It utilizes an optimization approach called **Variational Hybrid Quantum Annealing**, combining a hybrid quantum annealer (D-Wave) with a classical optimizer (SPSA).

---

## 📌 Table of Contents
- [About the Paper / Context](#-about-the-paper--context)
- [Compared QUBO Models](#-compared-qubo-models)
- [Project Architecture](#-project-architecture)
- [Installation & Prerequisites](#-installation--prerequisites)
- [API Keys Setup](#-api-keys-setup)
- [Usage](#-usage)
- [Methodology & Training (VHQAE)](#-methodology--training-vhqae)
- [Results & Evaluation](#-results--evaluation)
- [References & Citations](#-references--citations)

---

## ℹ️ About the Paper / Context

Predicting RNA secondary structure (including pseudoknots) under the minimum free energy (MFE) model is an **NP-hard** problem. 
This repository is based on the research presented in:

> **Title:** *A QUBO model of the RNA folding problem optimized by variational hybrid quantum annealing*  
> **Authors:** Tristan Zaborniak, Juan Giraldo, Hausi Müller, Hosna Jabbari, Ulrike Stege  
> **Published in:** 2022 IEEE International Conference on Quantum Computing and Engineering (QCE)  
> **DOI:** [10.1109/QCE53715.2022.00037](https://doi.org/10.1109/QCE53715.2022.00037)[cite: 1]

---

## 🧪 Compared QUBO Models

This project allows manipulating and comparing three distinct formulations:

1. **Model 1 (M1) - Baseline (Stem-level):** Maximizes base pairs and average stem length with heuristic penalties for pseudoknots.
2. **Model 2 (M2) - Stacked Quartets:** Based on length-2 sub-units and experimental nearest-neighbor stacking energies.
3. **Model 3 (M3) - Proposed Physical Model:** Incorporates stem thermodynamics (nearest-neighbor), a polymer physics-inspired penalty for pseudoknots ($P_{PK}'$), and a hairpin loop penalty.

---

## 📁 Project Architecture

```text
├── data/
│   ├── raw/                 # bpRNA-1m files (.ct / connectivity)
│   ├── train/               # Training dataset (70 structures)
│   └── test/                # Test dataset (40 structures)
├── src/
│   ├── preprocessing.py     # Stem & quartet extraction, base-pair matrices
│   ├── qubo_models.py       # Construction of Hamiltonians/QUBOs (M1, M2, M3)
│   ├── quantum_solver.py    # Interface with D-Wave Hybrid / Amazon Braket[cite: 1]
│   ├── optimizer.py         # Classical SPSA optimizer for parameters[cite: 1]
│   └── metrics.py           # Calculation of MCC (Matthews Correlation Coefficient) scores[cite: 1]
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   └── 02_results_analysis.ipynb
├── tests/                   # Unit tests
├── main_train.py            # VHQAE training loop[cite: 1]
├── main_test.py             # Test dataset evaluation script[cite: 1]
├── requirements.txt
└── README.md
