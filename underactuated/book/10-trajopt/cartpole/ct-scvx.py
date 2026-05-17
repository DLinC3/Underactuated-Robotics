"""Drake-style prox-linear minimum-time trajectory optimization.

This is a pydrake rewrite of the MATLAB Listing 4.2 minimum-time example.
The point is not to hide the optimization inside one long script.  The point is
Drake-style organization:

  1. A problem data object describes the world.
  2. A planner/optimizer computes a trajectory.
  3. A LeafSystem plant defines the continuous-time dynamics.
  4. A LeafSystem controller plays the optimized input trajectory.
  5. A Diagram wires plant + controller + logger.
  6. Simulator rolls out the closed-loop/feedforward system.

The dynamics are a 2D double integrator

    x = [px, py, vx, vy]
    u = [ax, ay]

with box-bounded acceleration, a speed bound, a circular obstacle, and an
unknown final time.  The final time is handled by time scaling: each interval
optimizes eta_k = [u_k; s_k], where s_k = dt / d tau.

Dependencies:
  - pydrake
  - numpy
  - matplotlib
  - a Drake-supported solver.  MOSEK is preferred for the conic speed bound.
    If the Lorentz cone binding is unavailable, the script falls back to a
    smooth nonlinear speed constraint, which typically requires SNOPT/IPOPT.

Run:
  python minimum_time_prox_linear_drake.py
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

import argparse
import math
import numpy as np

from pydrake.all import (
    BasicVector,
    DiagramBuilder,
    LeafSystem,
    LogVectorOutput,
    MathematicalProgram,
    MosekSolver,
    Simulator,
    Solve,
)


@dataclass
class MinimumTimeProblem:
    """All numerical data for the minimum-time example."""

    nx: int = 4
    nu: int = 2
    N: int = 10

    x0: np.ndarray = field(default_factory=lambda: np.zeros(4))
    xf: np.ndarray = field(default_factory=lambda: np.array([10.0, 10.0, 0.0, 0.0]))

    tmin: float = 5.0
    tmax: float = 20.0
    umax: np.ndarray = field(default_factory=lambda: np.ones(2))
    velmax: float = 3.0

    obst_cen: np.ndarray = field(default_factory=lambda: np.array([5.0, 4.5]))
    obst_rad: float = 3.0

    max_iter: int = 300
    tol: float = 5e-3
    gamma: float = 1e2
    eps_path: float = 1e-6

    fd_eps: float = 1e-5
    rk4_substeps: int = 20
    use_lorentz_speed_constraint: bool = True
    prefer_mosek: bool = True
    verbose: bool = True

    @property
    def rho(self) -> float:
        return 1.0 / self.gamma

    @property
    def tau(self) -> np.ndarray:
        return np.linspace(0.0, 1.0, self.N)

    @property
    def nxa(self) -> int:
        # augmented state: [physical state; accumulated path violation; time]
        return self.nx + 2

    @property
    def nua(self) -> int:
        # augmented input: [physical input; time-scale variable]
        return self.nu + 1

    @property
    def Ac(self) -> np.ndarray:
        return np.block(
            [[np.zeros((2, 2)), np.eye(2)], [np.zeros((2, 2)), np.zeros((2, 2))]]
        )

    @property
    def Bc(self) -> np.ndarray:
        return np.vstack([np.zeros((2, 2)), np.eye(2)])

    @property
    def umin(self) -> np.ndarray:
        return -self.umax

    @property
    def eta_lb(self) -> np.ndarray:
        return np.r_[self.umin, self.tmin]

    @property
    def eta_ub(self) -> np.ndarray:
        return np.r_[self.umax, self.tmax]


@dataclass
class ProxLinearIteration:
    iteration: int
    final_time: float
    slack: float
    step: float
    solver: str


@dataclass
class MinimumTimePlan:
    xi: np.ndarray              # shape: (nx + 2, N)
    eta: np.ndarray             # shape: (nu + 1, N - 1)
    physical_time: np.ndarray   # shape: (N,)
    iterations: List[ProxLinearIteration]

    @property
    def x_nodes(self) -> np.ndarray:
        return self.xi[:4, :]

    @property
    def u_intervals(self) -> np.ndarray:
        return self.eta[:2, :]

    @property
    def final_time(self) -> float:
        return float(self.physical_time[-1])


class ProxLinearMinimumTimePlanner:
    """Successive convexification / prox-linear planner.

    The MATLAB script is essentially this class written inline:

        current nominal trajectory
        -> linearize nonlinear flow map on every interval
        -> solve a convex prox-linear subproblem
        -> repeat

    Here the optimizer has a single public method, Solve(), because from the
    system point of view it is a trajectory source: it produces the plan that
    the controller will later play back.
    """

    def __init__(self, problem: MinimumTimeProblem):
        self.p = problem

    def Solve(self) -> MinimumTimePlan:
        xi_tilde, eta_tilde = self._initial_guess()
        iterations: List[ProxLinearIteration] = []

        for it in range(1, self.p.max_iter + 1):
            A, B, c = self._linearize_all_intervals(xi_tilde, eta_tilde)
            xi_opt, eta_opt, q_opt, r_opt, solver_name = self._solve_convex_subproblem(
                xi_tilde, eta_tilde, A, B, c
            )

            slack = float(np.sum(q_opt + r_opt))
            step = float(
                np.sqrt(np.sum((xi_opt - xi_tilde) ** 2) + np.sum((eta_opt - eta_tilde) ** 2))
            )
            time_grid = self._physical_time_from_eta(eta_opt)
            record = ProxLinearIteration(
                iteration=it,
                final_time=float(time_grid[-1]),
                slack=slack,
                step=step,
                solver=solver_name,
            )
            iterations.append(record)

            if self.p.verbose:
                print(
                    f"iter {it:03d} | T = {record.final_time:8.4f} | "
                    f"slack = {slack:9.3e} | step = {step:9.3e} | {solver_name}"
                )

            if step <= self.p.tol and slack <= self.p.tol:
                break

            xi_tilde = xi_opt
            eta_tilde = eta_opt

        return MinimumTimePlan(
            xi=xi_opt,
            eta=eta_opt,
            physical_time=self._physical_time_from_eta(eta_opt),
            iterations=iterations,
        )

    def _initial_guess(self) -> Tuple[np.ndarray, np.ndarray]:
        """A deliberately simple, smooth nominal trajectory.

        A straight line would pass directly through the obstacle.  We initialize
        with a high arc so that the first linearization is at least geometrically
        meaningful.  This is only a nominal trajectory; the optimizer is free to
        change it.
        """
        tau = self.p.tau
        T_guess = 10.0
        bump_height = 5.0

        xi = np.zeros((self.p.nxa, self.p.N))
        eta = np.zeros((self.p.nua, self.p.N - 1))

        for k, a in enumerate(tau):
            # cubic interpolation with zero endpoint derivative
            sigma = 3.0 * a**2 - 2.0 * a**3
            dsigma = 6.0 * a - 6.0 * a**2

            pos = (1.0 - sigma) * self.p.x0[:2] + sigma * self.p.xf[:2]
            pos[1] += bump_height * math.sin(math.pi * a) ** 2

            dpos_dtau = dsigma * (self.p.xf[:2] - self.p.x0[:2])
            dpos_dtau[1] += bump_height * math.pi * math.sin(2.0 * math.pi * a)
            vel = dpos_dtau / T_guess

            xi[:4, k] = np.r_[pos, vel]
            xi[4, k] = 0.0
            xi[5, k] = T_guess * a

        eta[:2, :] = 0.0
        eta[2, :] = T_guess
        return xi, eta

    def _linearize_all_intervals(
        self, xi_tilde: np.ndarray, eta_tilde: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        A = np.zeros((self.p.nxa, self.p.nxa, self.p.N - 1))
        B = np.zeros((self.p.nxa, self.p.nua, self.p.N - 1))
        c = np.zeros((self.p.nxa, self.p.N - 1))

        tau = self.p.tau
        for k in range(self.p.N - 1):
            A[:, :, k], B[:, :, k], c[:, k] = self._linearize_flow_map(
                xi_tilde[:, k], eta_tilde[:, k], tau[k], tau[k + 1]
            )
        return A, B, c

    def _linearize_flow_map(
        self, xbar: np.ndarray, eta: np.ndarray, tau0: float, tau1: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Linearize xbar_next = F_k(xbar, eta) by central finite differences.

        The MATLAB listing integrates variational equations to get the same
        objects.  This version uses finite differences for readability.  The
        subproblem still sees the same affine model

            xbar_{k+1} = A_k xbar_k + B_k eta_k + c_k.
        """
        F0 = self._flow_map_rk4(xbar, eta, tau0, tau1)
        A = np.zeros((self.p.nxa, self.p.nxa))
        B = np.zeros((self.p.nxa, self.p.nua))
        eps = self.p.fd_eps

        for i in range(self.p.nxa):
            d = np.zeros(self.p.nxa)
            d[i] = eps
            A[:, i] = (
                self._flow_map_rk4(xbar + d, eta, tau0, tau1)
                - self._flow_map_rk4(xbar - d, eta, tau0, tau1)
            ) / (2.0 * eps)

        for j in range(self.p.nua):
            d = np.zeros(self.p.nua)
            d[j] = eps
            B[:, j] = (
                self._flow_map_rk4(xbar, eta + d, tau0, tau1)
                - self._flow_map_rk4(xbar, eta - d, tau0, tau1)
            ) / (2.0 * eps)

        c = F0 - A @ xbar - B @ eta
        return A, B, c

    def _flow_map_rk4(self, xbar0: np.ndarray, eta: np.ndarray, tau0: float, tau1: float) -> np.ndarray:
        h = (tau1 - tau0) / self.p.rk4_substeps
        xbar = np.array(xbar0, dtype=float).copy()
        for _ in range(self.p.rk4_substeps):
            k1 = self._augmented_dynamics(xbar, eta)
            k2 = self._augmented_dynamics(xbar + 0.5 * h * k1, eta)
            k3 = self._augmented_dynamics(xbar + 0.5 * h * k2, eta)
            k4 = self._augmented_dynamics(xbar + h * k3, eta)
            xbar = xbar + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        return xbar

    def _augmented_dynamics(self, xbar: np.ndarray, eta: np.ndarray) -> np.ndarray:
        x = xbar[:4]
        pos = x[:2]
        vel = x[2:4]
        u = eta[:2]
        s = eta[2]

        speed_violation = np.dot(vel, vel) / self.p.velmax**2 - 1.0
        obstacle_violation = 1.0 - np.dot(pos - self.p.obst_cen, pos - self.p.obst_cen) / self.p.obst_rad**2
        path_integrand = 0.5 * max(speed_violation, 0.0) ** 2 + 0.5 * max(obstacle_violation, 0.0) ** 2

        xdot = self.p.Ac @ x + self.p.Bc @ u
        return s * np.r_[xdot, path_integrand, 1.0]

    def _solve_convex_subproblem(
        self,
        xi_tilde: np.ndarray,
        eta_tilde: np.ndarray,
        A: np.ndarray,
        B: np.ndarray,
        c: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, str]:
        prog = MathematicalProgram()
        xi = prog.NewContinuousVariables(self.p.nxa, self.p.N, "xi")
        eta = prog.NewContinuousVariables(self.p.nua, self.p.N - 1, "eta")
        q = prog.NewContinuousVariables(self.p.nxa, self.p.N - 1, "q")
        r = prog.NewContinuousVariables(self.p.nxa, self.p.N - 1, "r")

        # Boundary conditions.
        prog.AddBoundingBoxConstraint(np.r_[self.p.x0, 0.0, 0.0], np.r_[self.p.x0, 0.0, 0.0], xi[:, 0])
        prog.AddBoundingBoxConstraint(self.p.xf, self.p.xf, xi[:4, self.p.N - 1])

        # Nonnegative exact-penalty slacks.
        prog.AddBoundingBoxConstraint(0.0, np.inf, q.flatten())
        prog.AddBoundingBoxConstraint(0.0, np.inf, r.flatten())

        for k in range(self.p.N - 1):
            # Linearized dynamics with relaxation slacks:
            # xi_{k+1} = A_k xi_k + B_k eta_k + c_k + q_k - r_k.
            for i in range(self.p.nxa):
                expr = xi[i, k + 1] - c[i, k] - q[i, k] + r[i, k]
                for j in range(self.p.nxa):
                    expr -= A[i, j, k] * xi[j, k]
                for j in range(self.p.nua):
                    expr -= B[i, j, k] * eta[j, k]
                prog.AddLinearConstraint(expr, 0.0, 0.0)

            prog.AddBoundingBoxConstraint(self.p.eta_lb, self.p.eta_ub, eta[:, k])
            prog.AddBoundingBoxConstraint(-np.inf, self.p.eps_path, np.array([xi[4, k + 1]]))
            self._add_speed_constraint(prog, xi[2:4, k + 1])

        # Objective: final time + exact penalty + proximal term.
        prog.AddLinearCost(xi[5, self.p.N - 1])
        for k in range(self.p.N - 1):
            for i in range(self.p.nxa):
                prog.AddLinearCost(self.p.gamma * q[i, k])
                prog.AddLinearCost(self.p.gamma * r[i, k])

        prox = 0.0
        for k in range(self.p.N):
            for i in range(self.p.nxa):
                prox += (xi[i, k] - xi_tilde[i, k]) ** 2
        for k in range(self.p.N - 1):
            for i in range(self.p.nua):
                prox += (eta[i, k] - eta_tilde[i, k]) ** 2
        prog.AddQuadraticCost(0.5 / self.p.rho * prox)

        prog.SetInitialGuess(xi, xi_tilde)
        prog.SetInitialGuess(eta, eta_tilde)
        prog.SetInitialGuess(q, np.zeros_like(q, dtype=float))
        prog.SetInitialGuess(r, np.zeros_like(r, dtype=float))

        result, solver_name = self._solve_program(prog)
        if not result.is_success():
            raise RuntimeError(
                "The prox-linear subproblem did not solve successfully. "
                f"Solver: {solver_name}. Result: {result.get_solution_result()}"
            )

        return (
            result.GetSolution(xi),
            result.GetSolution(eta),
            result.GetSolution(q),
            result.GetSolution(r),
            solver_name,
        )

    def _add_speed_constraint(self, prog: MathematicalProgram, vel_vars: np.ndarray) -> None:
        """Add ||v||_2 <= v_max.

        The clean convex form is a Lorentz cone:

            [v_max, vx, vy] in Q_3.

        Drake exposes this when a conic solver path is available.  The fallback
        is a nonlinear scalar constraint; it is still the same feasible set, but
        it will not be treated as a conic constraint by generic solvers.
        """
        if self.p.use_lorentz_speed_constraint and hasattr(prog, "AddLorentzConeConstraint"):
            A = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
            b = np.array([self.p.velmax, 0.0, 0.0])
            try:
                prog.AddLorentzConeConstraint(A, b, vel_vars)
                return
            except TypeError:
                # Some older Drake wheels have a slightly different binding.
                pass

        def speed_margin(v: np.ndarray) -> np.ndarray:
            return np.array([self.p.velmax**2 - float(np.dot(v, v))])

        prog.AddConstraint(speed_margin, np.array([0.0]), np.array([np.inf]), vel_vars)

    def _solve_program(self, prog: MathematicalProgram):
        if self.p.prefer_mosek:
            try:
                mosek = MosekSolver()
                if mosek.available() and mosek.enabled():
                    return mosek.Solve(prog), "MOSEK"
            except Exception:
                pass
        result = Solve(prog)
        solver_id = result.get_solver_id()
        solver_name = solver_id.name() if solver_id is not None else "Drake default solver"
        return result, solver_name

    def _physical_time_from_eta(self, eta: np.ndarray) -> np.ndarray:
        ts = np.zeros(self.p.N)
        tau = self.p.tau
        for k in range(self.p.N - 1):
            ts[k + 1] = ts[k] + eta[2, k] * (tau[k + 1] - tau[k])
        return ts


class DoubleIntegrator2D(LeafSystem):
    """A tiny Drake System: xdot = Ac x + Bc u."""

    def __init__(self):
        super().__init__()
        self._Ac = np.block(
            [[np.zeros((2, 2)), np.eye(2)], [np.zeros((2, 2)), np.zeros((2, 2))]]
        )
        self._Bc = np.vstack([np.zeros((2, 2)), np.eye(2)])

        self.DeclareVectorInputPort("u", BasicVector(2))
        self.DeclareContinuousState(4)
        self.DeclareVectorOutputPort("x", BasicVector(4), self._copy_state_out)

    def _copy_state_out(self, context, output):
        output.SetFromVector(context.get_continuous_state_vector().CopyToVector())

    def DoCalcTimeDerivatives(self, context, derivatives):
        x = context.get_continuous_state_vector().CopyToVector()
        u = self.get_input_port(0).Eval(context)
        derivatives.get_mutable_vector().SetFromVector(self._Ac @ x + self._Bc @ u)


class PiecewiseConstantFeedforwardController(LeafSystem):
    """A controller System that plays the optimized input trajectory.

    In a more serious robot problem this block is where we would put feedback,
    tracking, saturation, or a state machine.  Here it is intentionally simple:
    the optimized trajectory is the policy.
    """

    def __init__(self, knot_times: np.ndarray, controls: np.ndarray):
        super().__init__()
        assert controls.shape[0] == 2
        assert controls.shape[1] == len(knot_times) - 1
        self._times = np.asarray(knot_times, dtype=float)
        self._controls = np.asarray(controls, dtype=float)
        self.DeclareVectorOutputPort("u", BasicVector(2), self._calc_control)

    def _calc_control(self, context, output):
        t = context.get_time()
        k = int(np.searchsorted(self._times, t, side="right") - 1)
        k = int(np.clip(k, 0, self._controls.shape[1] - 1))
        output.SetFromVector(self._controls[:, k])


def make_diagram(plan: MinimumTimePlan):
    """Build a tiny Drake Diagram: controller -> plant -> logger."""
    builder = DiagramBuilder()

    plant = builder.AddSystem(DoubleIntegrator2D())
    controller = builder.AddSystem(
        PiecewiseConstantFeedforwardController(plan.physical_time, plan.u_intervals)
    )

    builder.Connect(controller.get_output_port(0), plant.get_input_port(0))
    logger = LogVectorOutput(plant.get_output_port(0), builder)

    diagram = builder.Build()
    diagram.set_name("minimum_time_double_integrator_diagram")
    return diagram, plant, logger


def simulate_plan(problem: MinimumTimeProblem, plan: MinimumTimePlan):
    diagram, plant, logger = make_diagram(plan)
    context = diagram.CreateDefaultContext()

    plant_context = plant.GetMyMutableContextFromRoot(context)
    plant_context.SetContinuousState(problem.x0)

    simulator = Simulator(diagram, context)
    simulator.Initialize()
    simulator.AdvanceTo(plan.final_time)

    log = logger.FindLog(context)
    return log.sample_times(), log.data()


def plot_plan_and_rollout(problem: MinimumTimeProblem, plan: MinimumTimePlan, t: np.ndarray, x: np.ndarray):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 7))
    theta = np.linspace(0.0, 2.0 * np.pi, 300)
    obs_x = problem.obst_cen[0] + problem.obst_rad * np.cos(theta)
    obs_y = problem.obst_cen[1] + problem.obst_rad * np.sin(theta)
    ax.fill(obs_x, obs_y, alpha=0.3, label="obstacle")
    ax.plot(plan.x_nodes[0, :], plan.x_nodes[1, :], "o--", label="optimized nodes")
    ax.plot(x[0, :], x[1, :], linewidth=2, label="Drake rollout")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)
    ax.set_xlabel("x position")
    ax.set_ylabel("y position")
    ax.set_title(f"Minimum-time trajectory, T = {plan.final_time:.3f} s")
    ax.legend()

    speed = np.sqrt(np.sum(x[2:4, :] ** 2, axis=0))
    node_speed = np.sqrt(np.sum(plan.x_nodes[2:4, :] ** 2, axis=0))
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    ax2.plot(t, speed, linewidth=2, label="Drake rollout")
    ax2.plot(plan.physical_time, node_speed, "o--", label="optimized nodes")
    ax2.axhline(problem.velmax, linestyle="--", label="speed limit")
    ax2.grid(True)
    ax2.set_xlabel("time")
    ax2.set_ylabel("speed")
    ax2.set_title("Speed profile")
    ax2.legend()
    return fig, fig2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-show", action="store_true", help="Do not call matplotlib show().")
    parser.add_argument("--save-plots", action="store_true", help="Save figures as PNG files.")
    parser.add_argument("--quiet", action="store_true", help="Suppress iteration printout.")
    args = parser.parse_args()

    problem = MinimumTimeProblem(verbose=not args.quiet)
    planner = ProxLinearMinimumTimePlanner(problem)
    plan = planner.Solve()

    print("\nPlanning complete")
    print(f"  final time: {plan.final_time:.6f}")
    print(f"  iterations: {len(plan.iterations)}")
    print(f"  final slack: {plan.iterations[-1].slack:.3e}")
    print(f"  final step:  {plan.iterations[-1].step:.3e}")

    t, x = simulate_plan(problem, plan)
    figs = plot_plan_and_rollout(problem, plan, t, x)

    if args.save_plots:
        figs[0].savefig("minimum_time_trajectory.png", dpi=200, bbox_inches="tight")
        figs[1].savefig("minimum_time_speed.png", dpi=200, bbox_inches="tight")

    if not args.no_show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()
