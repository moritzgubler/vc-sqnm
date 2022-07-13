# vc-sqnm
Implementation of the vc-sqnm and sqnm optimization algorithms in C++

The stabilized quasi Newton method (SQNM) is a fast and reliable optimization method that is well adapted to find local minima on the potential energy surface. When using the SQNM method please cite [the initial SQNM paper](https://aip.scitation.org/doi/10.1063/1.4905665).

This repository also contains an implementation of the variable cell shape SQNM method which optimizes also the lattice vectors of systems with periodic boundary conditions. When the vc-SQNM method is used please cite the [vc-SQNM](https://arxiv.org/abs/2206.07339) paper as well. 
