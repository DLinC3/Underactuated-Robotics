# barrierQP

A tiny Mehrotra-style predictor-corrector primal-dual interior-point QP solver. It rebuilds method used in solvers CVXGEN/qpSWIFT/HPIPM from scratch: affine prediction, adaptive centering/correction, positivity-preserving step control, strictly interior initialization, and residual-based safeguards. Useful for seeing the algorithmic logic behind CVXGEN/qpSWIFT/HPIPM QP solvers.

## Clone

To clone only the standalone solver:

```bash
git clone https://github.com/denglinc/barrierQP.git
```