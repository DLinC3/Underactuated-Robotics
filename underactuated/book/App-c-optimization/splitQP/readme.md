# splitQP

A tiny JAX solver for families of dense convex quadratic programs that share a
fixed quadratic and constraint matrix. It implements proximal ADMM, and its one
idea is factorization reuse: `Solver(P, A)` factors a single matrix at
construction and reuses it for every iteration and every member of the family.
This is a compact educational research implementation, not a production solver.

## Clone

To clone only the standalone solver:

```bash
git clone https://github.com/denglinc/splitQP.git
```
