# EE 521: Analysis of Power Systems and EE 523: Power System Stability and Control
## Fall 2022 | Noel Schulz and Spring 2023 | Mani V. Venkatasubramanian 
 <img src = "https://user-images.githubusercontent.com/24756405/237030072-5b15d383-9fec-4af8-bf5d-572b6db31e37.png" width = 35% height = 30%> <img src = "https://user-images.githubusercontent.com/24756405/237032599-14edd7bc-5b4c-4f0d-84a3-3ff7e9f963ec.png" width = 50% height = 100%>

### MATLAB implementations for the two courses at Washington State University, Pullman.

### Power Flow Algorithms added:
- Newton Raphson Power Flow `NRPF`
- `Decoupled NRPF`
- `Fast Decoupled NRPF`
- `Continuation Power Flow`

### Linear System Solving Algorithms added:
- `LU Factorization`

### Textbook solved examples added:
- `koth3`: A 3 bus system from Kothari and Nagrath's Modern Power System Analysis.
- `crow3`: The 3 bus system in Example 5.9 from Mariesa L Crow's Computational Methods for Electric Power Systems. 
- `ieee11`: Kundur's 2 Area 11 bus system as given in Example 12.6, Pg 813 of Power System Stability and Control by Prabha Kundur.
### Data Strucutures and Algorithms Sparsified:
- `YBus`
- Jacobian `J`
- Computation of Mismatches $[\Delta P ;\Delta Q]$.
- `sparmat` and `sparvec` can convert matrices and vectors in compressed format `(nrow, ncol, val)` or `(nIndex, val)` into the sparse format `[nnz, N]`. All data structures are tables.

## Stability and Control Scripts for the Kundur 4 Machine 2 Area System:
| Model Type      | Dynamic Initialization | Small Signal Stability Analysis   | Transient Stability Analysis |
| :---        |    :----:   |          :---: | :---: |
| Type 3 aka Classical Model      | ✅       | ✅   | ✅ |
| Type 2 with AVR and Governor   | ✅        | 🟨      | 🟨 |
| Type 1 with AVR and Governor | ✅ | ✅ | 🔴 | 

Legend:
| Symbol | Remark |
| :---: | :--- |
| ✅ | Implemented and performing as expected |
| 🟨 | Implemented but NOT performing as expected |
| 🔴 | NOT implemented |

### Yet to implement:
- Sparse LU Factorization
- [OPTIONAL] Bus Changing (PV to PQ)
- [OPTIONAL] DC Power Flow

### Caveats:
- Currently the NRPF algorithm does NOT support bus type conversion. This obviously affects convergence for bigger bus systems, but fortunately does not seem to affect the `ieee14` and `ieee30` bus systems.
- Currently it converges for the `ieee14` and `ieee30` bus systems, has trouble with `ieee57` bus system for a couple of buses and blows up for the `ieee118` bus system.
- $N$ bus systems with individual bus numbers $i$ outside the range of natural numbers from $[1, N]$  are currently NOT supported. The `ieee300` bus system is one such system.

<img src = "https://user-images.githubusercontent.com/24756405/237028282-16ffc864-98f8-4c6f-a663-dfb04d191623.png" width = 15% height = 15%>
