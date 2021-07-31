# d3_qpgf_d1
This is the calculation program for quasi-periodic Green's function of the Helmholtz equations. The quasi-periodicity is 1-dimension ( x axis ), Green's function is 3-dimensions. This program is used my original method.

## Definitions
- quasi-periodic Green's function

  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r})=\sum_{l=-\infty}^{\infty}G(\mathbf{r}+l\mathbf{d})\exp(il\mathbf{k}\cdot\mathbf{d})">　　  
  <img src="https://latex.codecogs.com/gif.latex?G(\mathbf{r})=\frac{\exp(ik|\mathbf{r}|)}{4\pi|\mathbf{r}|}">  
  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r})"> is quasi-periodic Green's function  
  <img src="https://latex.codecogs.com/gif.latex?G(\mathbf{r})"> is Green's function of the 3-dimensional Helmholtz equation  
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{k}"> is wave number vector, 
  <img src="https://latex.codecogs.com/gif.latex?|\mathbf{k}|=k">  
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{d}"> is lattice vector, where
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{d}=(d,0,0)">

- quasi-peridic condition

  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r}+l\mathbf{d})=\exp(-il\mathbf{k}\cdot\mathbf{d})\,^q\!G(\mathbf{r})">

## Usage of example code
1. type 'make' command to compile
2. type './example.out' to run

Please see src/d3_qpgf_d1.h for detail of functions, src/example.c for detail of function usages.

## References
1. Capolino, Filippo, Donald R. Wilton, and William A. Johnson. "Efficient computation of the 2-D Green's function for 1-D periodic structures using the Ewald method." IEEE Transactions on Antennas and Propagation 53.9 (2005): 2977-2984.  
2. Beylkin, Gregory, Christopher Kurcz, and Lucas Monzón. "Fast algorithms for Helmholtz Green's functions." Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 464.2100 (2008): 3301-3326.
3. Abramowitz, Milton, Irene A. Stegun, and Robert H. Romer. "Handbook of mathematical functions with formulas, graphs, and mathematical tables." (1988): 958-958.
4. Faddeeva Package. http://ab-initio.mit.edu/Faddeeva
