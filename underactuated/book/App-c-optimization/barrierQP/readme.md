# barrierQP

A tiny JAX solver for dense, strictly convex quadratic programs with equality and
inequality constraints. It implements Mehrotra's primal-dual predictor-corrector
interior-point method, and its one idea is factorization reuse: each iteration
factors a single KKT matrix and reuses it for both Newton directions. This is a
compact educational research implementation, not a production solver.

## Clone

To clone only the standalone solver:

```bash
git clone https://github.com/denglinc/barrierQP.git
```
