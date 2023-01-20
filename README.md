# EE 521: Analysis of Power Systems and EE 523: Power System Stability and Control
## Fall 2022 | Noel Schulz and Spring 2023 | Mani V. Venkatasubramanian
### MATLAB implementations for the two courses at Washington State University, Pullman.

- Make sure to set MATLAB's working directory as the location of ee521_and_ee523_main.mlx
- Currently the NRPF algorithm does NOT support bus type conversion. This obviously affects convergence for bigger bus systems, but fortunately does not seem to affect the `ieee14` and `ieee30` bus systems.
- Currently it converges for the `ieee14` and `ieee30` bus systems, has trouble with `ieee57` bus system for a couple of buses and blows up for the `ieee118` bus system.
- $N$ bus systems with individual bus numbers $i$ outside the range of natural numbers from $[1, N]$  are currently NOT supported. The `ieee300` bus system is one such system.

### Power Flow Algorithms added:
- Newton Raphson Power Flow `NRPF`
- `Decoupled NRPF`
- `Fast Decoupled NRPF`

### Linear System Solving Algorithms added:
- `LU Factorization`

### Textbook solved examples added:
- `koth3`: A 3 bus system from Kothari and Nagrath's Modern Power System Analysis.
- `crow3`: The 3 bus system in Example 5.9 from Mariesa L Crow's Computational Methods for Electric Power Systems. 

### Data Strucutures and Algorithms Sparsified:
- `YBus`
- Jacobian `J`
- Computation of Mismatches $[\Delta P ;\Delta Q]$.
- `sparmat` and `sparvec` can convert matrices and vectors in compressed format `(nrow, ncol, val)` or `(nIndex, val)` into the sparse format `[nnz, N]`. All data structures are tables.

### Yet to implement:
- Sparse LU Factorization
- [OPTIONAL] Bus Changing (PV to PQ)
- [OPTIONAL] DC Power Flow
