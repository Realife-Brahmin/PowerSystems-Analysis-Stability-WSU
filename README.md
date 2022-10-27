# EE 521: Power System Analysis
## Fall 2022 | Noel Schulz
### MATLAB implementations for the Power System Analysis course (EE 521) at WSU.

- Make sure to set MATLAB's working directory as the location of ee521_powerFlow.mlx
- Currently the NRPF algorithm does NOT support bus type conversion. This obviously affects convergence for bigger bus systems, but fortunately does not seem to affect the `ieee14` and `ieee30` bus systems.
- Currently it converges for the `ieee14` and `ieee30` bus systems, has trouble with `ieee57` bus system for a couple of buses and blows up for the `ieee118` bus system. I have added another small 3 bus system which I am calling `koth3` (from Kothari and Nagrath's Modern Power System Analysis) for sanity checks, which converges. 
- $N$ bus systems with individual bus numbers $i$ outside the range of natural numbers from $[1, N]$  are currently NOT supported. The `ieee300` bus system is one such system.
