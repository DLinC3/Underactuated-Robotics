# splitQP

A tiny OSQP-style QP solver. It rebuilds the ADMM / operator-splitting core from scratch: augmented-Lagrangian splitting, sparse KKT solves, over-relaxation, factorization caching, adaptive rho updates, and residual-based stopping. Small enough to read, but faithful enough to show OSQP-style iterations.

## Clone

To clone only the standalone solver:

```bash
git clone https://github.com/denglinc/splitQP.git
```