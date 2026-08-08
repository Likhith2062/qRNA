# 🧬 RNA Secondary Structure Prediction using Quantum Annealing

<div align="center">

### A Modular Research Implementation of RNA Secondary Structure Prediction using
### **Quadratic Unconstrained Binary Optimization (QUBO)** and **Quantum Annealing**

*Built using ViennaRNA, D-Wave Ocean SDK and Python*

![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)
![ViennaRNA](https://img.shields.io/badge/ViennaRNA-2.x-green.svg)
![D-Wave Ocean](https://img.shields.io/badge/D--Wave-Ocean-purple.svg)
![Optimization](https://img.shields.io/badge/Optimization-QUBO-orange.svg)
![Platform](https://img.shields.io/badge/Platform-qBraid%20%7C%20Local-black.svg)

</div>

---

> [!IMPORTANT]
> This repository contains a complete implementation of an RNA secondary structure prediction framework based on **Quadratic Unconstrained Binary Optimization (QUBO)**. The project reformulates RNA folding as a combinatorial optimization problem that can be solved using **quantum-inspired annealing** or executed on compatible **quantum annealing hardware** with minimal modifications.

---

# 📖 Abstract

Predicting the secondary structure of RNA is one of the fundamental problems in computational biology. Since RNA molecules fold into energetically favorable configurations through intramolecular base pairing, determining the correct secondary structure enables a deeper understanding of numerous biological processes including gene regulation, protein synthesis, catalysis, viral replication and RNA-based therapeutics.

Unlike traditional dynamic programming approaches that explicitly search the folding space, this project reformulates RNA folding as a **Quadratic Unconstrained Binary Optimization (QUBO)** problem. Candidate stems are first generated geometrically, thermodynamic information is obtained from the ViennaRNA Package, and the complete folding problem is encoded into a Hamiltonian suitable for optimization through quantum annealing.

The implementation follows the formulation proposed by **Zaborniak et al.**, while extending the Hamiltonian with a polymer entropy based pseudoknot penalty inspired by the **ShapeKnots** model proposed by **Hajdin et al.** The resulting optimization problem is solved using D-Wave Ocean's simulated annealing backend, while maintaining compatibility with actual quantum annealers through an abstract solver interface.

This repository serves both as a research implementation and as an educational reference illustrating how modern RNA thermodynamic models can be translated into optimization problems suitable for quantum computing.

---

# ✨ Project Highlights

✔ Complete implementation of the Model-3 QUBO formulation proposed by Zaborniak et al.

✔ ViennaRNA integration for nearest-neighbour thermodynamic calculations

✔ Automatic geometric stem generation

✔ Duplicate-free candidate stem enumeration

✔ Polymer entropy based pseudoknot modelling

✔ Quantum-inspired optimization using D-Wave Ocean SDK

✔ Modular architecture supporting future quantum hardware execution

✔ Automatic decoding into dot-bracket notation

✔ Evaluation using Matthews Correlation Coefficient (MCC)

✔ Comprehensive unit tests for every major module

---

# 📌 Table of Contents

- [Project Motivation](#-project-motivation)
- [Repository Structure](#-repository-structure)
- [Scientific Background](#-scientific-background)
- [Mathematical Formulation](#-mathematical-formulation)
- [Software Architecture](#-software-architecture)
- [Installation](#-installation)
- [Execution](#-execution)
- [Classical Benchmark Results](#-classical-benchmark-results)
- [Quantum / Quantum-inspired Implementation](#-quantum--quantum-inspired-implementation)
- [Results and Analysis](#-results-and-analysis)
- [Scaling and Quantum Resource Analysis](#-scaling-and-quantum-resource-analysis)
- [Assumptions](#-assumptions)
- [Limitations](#-limitations)
- [Future Work](#-future-work)
- [References](#-references)

---

# 🎯 Project Motivation

RNA secondary structure prediction is widely regarded as an NP-hard optimization problem due to the exponentially increasing number of possible folding configurations as sequence length grows.

Traditional RNA folding algorithms generally rely on dynamic programming under the Minimum Free Energy (MFE) framework. Although these approaches are computationally efficient for many practical problems, incorporating complex structural motifs such as pseudoknots substantially increases computational complexity.

Recent advances in quantum optimization provide an alternative perspective by reformulating RNA folding as a combinatorial optimization problem. Rather than explicitly constructing every possible secondary structure, the problem is represented as a binary optimization task where each candidate stem corresponds to a binary decision variable. The globally optimal RNA structure is then obtained by minimizing a Hamiltonian encoding thermodynamic stability and structural constraints.

The primary goal of this project is to investigate this optimization-based formulation while producing a modular software implementation suitable for experimentation, benchmarking and future execution on quantum annealing hardware.

---

# 📁 Repository Structure

```text
Project
│
├── main.py
│   End-to-end RNA secondary structure prediction pipeline.
│
├── vienna.py
│   Thin wrapper around ViennaRNA providing thermodynamic
│   lookup functions and structure evaluation utilities.
│
├── stem_generator.py
│   Geometric candidate stem generation.
│
├── qubo_builder.py
│   Construction of the complete RNA folding Hamiltonian.
│
├── quantum_solver.py
│   QUBO optimization using D-Wave Ocean SDK.
│
├── decoder.py
│   Converts binary optimization variables into
│   RNA secondary structures.
│
├── metrics.py
│   Evaluation metrics including Matthews Correlation
│   Coefficient (MCC).
│
├── tests/
│   Comprehensive unit tests.
│
├── experiments/
│   Validation and exploratory scripts used during development.
│
└── README.md
```

---

# 🏗 Overall Pipeline

```text
                     RNA Sequence
                          │
                          ▼
                  ViennaRNA Wrapper
                          │
                          ▼
                Candidate Stem Generator
                          │
                          ▼
                 Candidate Stem Library
                          │
                          ▼
                  QUBO Hamiltonian Builder
                          │
                          ▼
            Quadratic Optimization Problem
                          │
                          ▼
         Quantum / Quantum-inspired Annealer
                          │
                          ▼
                Optimal Binary Solution
                          │
                          ▼
                    Structure Decoder
                 ┌────────┴────────┐
                 ▼                 ▼
          Dot-Bracket       Adjacency Matrix
                 │                 │
                 └────────┬────────┘
                          ▼
            ViennaRNA Benchmark Comparison
                          │
                          ▼
          Matthews Correlation Coefficient
```

---

> [!NOTE]
> The implementation intentionally separates **thermodynamic calculations**, **candidate stem generation**, **QUBO construction**, **optimization**, **structure decoding**, and **evaluation** into independent modules. This modular design simplifies testing, maintenance, future extensions and replacement of individual components (e.g., executing on real D-Wave quantum hardware instead of a simulated annealer).