# Project Title: Numerical Solution of 2D Steady-State Diffusion Equation

## Overview

This project aims to numerically solve the 2D steady-state diffusion equation using various numerical methods, including:

1. **Gauss Elimination**
2. **Gauss-Seidel Iteration**
3. **Line-by-Line (Row Sweep)**
4. **Alternating Direction Implicit (ADI)**

The project includes MATLAB implementations of these methods, along with code to visualize the solutions and analyze their performance.

## Implementation Details

**main.m:**

* Initializes the domain and boundary conditions.
* Calls the respective solvers (Gauss elimination, Gauss-Seidel, row sweep, column sweep, and ADI).
* Visualizes the solutions using contour plots.

**gauss_elimination.m:**

* Implements the Gauss elimination method to solve the system of linear equations.
* Converts the coefficient matrix to row-echelon form.
* Performs back-substitution to obtain the solution.

**gauss_seidel.m:**

* Implements the Gauss-Seidel iterative method.
* Iteratively updates the solution values using the latest available values.

**row_sweep.m:**

* Implements the line-by-line (row sweep) method.
* Solves for the unknowns in each row using the Thomas algorithm (tridiagonal matrix algorithm).

**column_sweep.m:**

* Implements the line-by-line (column sweep) method.
* Solves for the unknowns in each column using the Thomas algorithm.

**ADI.m:**

* Implements the Alternating Direction Implicit (ADI) method.
* Alternates between row-wise and column-wise sweeps to solve the implicit equations.

## Results and Analysis

* **Contour plots:** Visualize the numerical solutions obtained from each method.
* **Convergence analysis:** Analyze the convergence rate and stability of each method.
* **Computational efficiency:** Compare the CPU time and memory usage of different methods.
