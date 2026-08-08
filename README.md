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

---

# 🧬 Scientific Background

## 🧪 What is RNA?

Ribonucleic Acid (RNA) is one of the fundamental biomolecules responsible for storing, transmitting and regulating genetic information inside living organisms. Unlike DNA, which typically exists as a stable double-stranded helix, RNA molecules are generally single stranded and capable of folding into highly complex three-dimensional conformations.

The biological function of an RNA molecule is determined not only by its nucleotide sequence but also by the secondary and tertiary structures formed through intramolecular base pairing.

RNA secondary structures play critical roles in numerous biological processes including:

- Protein synthesis (mRNA, tRNA and rRNA)
- Gene regulation
- RNA interference
- Ribozymes and catalytic RNAs
- Viral genome organization
- Drug discovery and RNA therapeutics

Accurate prediction of RNA secondary structure therefore remains one of the central problems in computational biology.

---

## 🔬 RNA Secondary Structure

RNA secondary structure describes the collection of hydrogen bonds formed between complementary nucleotides within the same RNA molecule.

The canonical Watson-Crick base pairs are

- Guanine — Cytosine (G-C)
- Adenine — Uracil (A-U)

along with the biologically important wobble pair

- Guanine — Uracil (G-U)

These interactions produce structural motifs such as

- Hairpin loops
- Internal loops
- Bulges
- Multibranch loops
- Helices
- Pseudoknots

A complete RNA secondary structure consists of a combination of these motifs arranged to minimize the overall free energy of the molecule.

---

## 🌡 Thermodynamic Folding

RNA folding is generally modeled using the **Minimum Free Energy (MFE)** principle.

The free energy of a secondary structure is determined using experimentally measured nearest-neighbour thermodynamic parameters, commonly referred to as the Turner Energy Rules.

The total free energy of an RNA structure is composed of several contributions:

- Base stacking interactions
- Hairpin loop penalties
- Internal loop penalties
- Bulge penalties
- Multiloop penalties
- Dangling end corrections

Among these, nearest-neighbour stacking interactions contribute the majority of structural stability.

This project uses the **ViennaRNA Package** as the thermodynamic engine for evaluating these energetic contributions.

---

## 📦 Candidate Stem Representation

Rather than searching directly over every possible RNA secondary structure, this implementation first generates every geometrically valid candidate stem.

Each candidate stem is defined by

- Starting position on the 5' strand
- Starting position on the 3' strand
- Stem length
- Total nearest-neighbour stacking energy
- Hairpin loop penalty

Each stem becomes a binary optimization variable within the QUBO Hamiltonian.

If

```
qi = 1
```

the stem is included in the predicted structure.

Otherwise

```
qi = 0
```

the stem is excluded.

This transformation converts RNA folding into a binary optimization problem.

---

# ⚛ Why QUBO?

Quadratic Unconstrained Binary Optimization (QUBO) is one of the most widely studied optimization models in quantum computing.

A generic QUBO problem is written as

\[
\min_x \; x^TQx
\]

where

- **Q** is a symmetric matrix of optimization coefficients.
- **x** is a binary vector.

Every binary variable can represent a decision within the optimization problem.

In this implementation,

```
qi
```

represents the decision

> "Should candidate stem *i* be included in the final RNA structure?"

The objective of the optimization process is therefore to identify the subset of stems that minimizes the Hamiltonian while satisfying all structural constraints.

---

## ⚙ Why Quantum Annealing?

Quantum annealing is an optimization technique designed specifically for solving QUBO and Ising Hamiltonian problems.

Instead of exhaustively searching every possible RNA folding configuration, the optimization algorithm searches the energy landscape defined by the Hamiltonian.

Advantages include:

- Native support for binary optimization problems.
- Natural mapping from Hamiltonians to hardware.
- Scalability to increasingly complex optimization formulations.
- Compatibility with D-Wave quantum annealers.

Although this project currently employs D-Wave Ocean's **Simulated Annealing Sampler**, the software architecture has been intentionally designed such that the optimization backend may be replaced with an actual quantum annealer with minimal modifications.

---

# 📚 Research Foundation

The mathematical formulation implemented throughout this repository is based primarily upon two research publications.

### Paper I

> **A QUBO Model of the RNA Folding Problem Optimized by Variational Hybrid Quantum Annealing**

This publication introduces three successive QUBO formulations for RNA folding.

The present implementation adopts **Model 3**, which incorporates

- nearest-neighbour stacking energies,
- hairpin loop penalties,
- overlap constraints,
- trainable Hamiltonian coefficients,
- and pseudoknot penalties.

---

### Paper II

> **Accurate SHAPE-directed RNA Secondary Structure Modeling, Including Pseudoknots**

The original QUBO formulation employs heuristic pseudoknot penalties.

To provide a more physically meaningful representation, this implementation instead adopts the polymer entropy model proposed by Hajdin et al.

The pseudoknot penalty therefore depends upon

- helix lengths,
- connecting loop lengths,
- polymer entropy,
- and empirically derived λ lookup tables,

rather than a constant heuristic penalty.

---

# 🧠 Design Philosophy

During development, a major design objective was strict separation between

1. Thermodynamic calculations
2. Candidate generation
3. Mathematical formulation
4. Optimization
5. Decoding
6. Evaluation

Consequently,

- **ViennaRNA** is responsible solely for thermodynamic calculations.
- **StemGenerator** performs only geometric candidate generation.
- **QUBOBuilder** constructs the Hamiltonian.
- **QuantumSolver** performs optimization.
- **Decoder** converts binary variables into RNA structures.
- **Metrics** evaluates prediction quality.

This separation significantly improves maintainability, unit testing and future extensibility.

---

> [!TIP]
> The implementation intentionally avoids using ViennaRNA to predict the final RNA structure. ViennaRNA is employed exclusively for thermodynamic parameter extraction and benchmarking. The final predicted secondary structure is obtained solely by minimizing the constructed QUBO Hamiltonian using a quantum-inspired optimization backend.