## Underactuated Robotics

This repository collects my worked notebooks, implementation notes, and supplementary code while studying Underactuated Robotics. Here are a few summaries:

---

<details>
<summary><strong>Drake Map</strong></summary>

Drake combines the Systems Framework, multibody modeling, visualization, simulation, and optimization into one coherent robotics workflow (Systems, MultibodyPlant, and MathematicalProgram.).

- [`underactuated/book/00-tutorials/index.ipynb`](underactuated/book/00-tutorials/index.ipynb)  
  A high-level index of the official Drake tutorials in this repository.

#### 1. Systems Framework

This section reviews `LeafSystem`, input/output ports, state, parameters, events, and `Context`.

- [`underactuated/book/00-tutorials/dynamical_systems.ipynb`](underactuated/book/00-tutorials/dynamical_systems.ipynb)  
  Introduces dynamical systems, simulation, `Context`, logging, and simple diagrams.

- [`underactuated/book/00-tutorials/authoring_leaf_systems.ipynb`](underactuated/book/00-tutorials/authoring_leaf_systems.ipynb)  
  The main reference for writing custom `LeafSystem`s with ports, state, parameters, events, and scalar conversion.

- [`underactuated/book/01-intro/exercises/drake_systems.ipynb`](underactuated/book/01-intro/exercises/drake_systems.ipynb)  
  A compact exercise for implementing a simple Drake system from scratch.

---

#### 2. Diagrams, Composition, and Contexts

This section reviews `DiagramBuilder`, `Diagram`, subsystem contexts, exported ports, and system composition.

- [`underactuated/book/00-tutorials/working_with_diagrams.ipynb`](underactuated/book/00-tutorials/working_with_diagrams.ipynb)  
  The cleanest tutorial for building diagrams, nesting diagrams, exporting ports, and accessing subsystem contexts.

- [`underactuated/book/08-lqr/exercises/drake_diagrams.ipynb`](underactuated/book/08-lqr/exercises/drake_diagrams.ipynb)  
  Uses a simple actuator and double-integrator example to practice diagram construction and LQR design.

- [`underactuated/book/02-pend/exercises/vibrating_pendulum.ipynb`](underactuated/book/02-pend/exercises/vibrating_pendulum.ipynb)  
  An integrated example connecting a multibody plant, controller, visualizer, simulator, and logger.

---

#### 3. Multibody Modeling and Simulation

This section reviews URDF/SDF parsing, `MultibodyPlant`, `SceneGraph`, geometry, and `Finalize()`.

- [`underactuated/book/00-tutorials/authoring_multibody_simulation.ipynb`](underactuated/book/00-tutorials/authoring_multibody_simulation.ipynb)  
  The main tutorial for loading robot models, building `MultibodyPlant + SceneGraph`, finalizing the plant, and running simulation.

- [`underactuated/book/03-acrobot/exercises/cartpoles_urdf.ipynb`](underactuated/book/03-acrobot/exercises/cartpoles_urdf.ipynb)  
  A hands-on URDF exercise for building, parsing, visualizing, and controlling cart-pole models.

- [`underactuated/book/00-tutorials/mathematical_program_multibody_plant.ipynb`](underactuated/book/00-tutorials/mathematical_program_multibody_plant.ipynb)  
  Builds an optimization problem around an IIWA `MultibodyPlant`, especially for IK-style problems, using custom evaluators compatible with both `float` and `AutoDiffXd`.

---

#### 4. Visualization, Animation, and Rendering

This section separates basic visualization, animation playback, interactive Meshcat controls, camera rendering, and lighting.

- [`underactuated/book/00-tutorials/rendering_multibody_plant.ipynb`](underactuated/book/00-tutorials/rendering_multibody_plant.ipynb)  
  The main reference for camera rendering, RGB-D images, label images, render engines, and `SceneGraphInspector`.

- [`underactuated/book/03-acrobot/cartpole.ipynb`](underactuated/book/03-acrobot/cartpole.ipynb)  
  The best notebook for interactive Meshcat usage, including sliders, buttons, gamepad teleoperation, and `MeshcatSliders`.

- [`underactuated/book/00-tutorials/pyplot_animation_multibody_plant.ipynb`](underactuated/book/00-tutorials/pyplot_animation_multibody_plant.ipynb)  
  Shows how to create planar PyPlot animations and notebook playback from a simulated multibody system.

- [`underactuated/book/00-tutorials/configuring_rendering_lighting.ipynb`](underactuated/book/00-tutorials/configuring_rendering_lighting.ipynb)  
  Covers advanced rendering configuration, lighting, and render-engine settings.

---

#### 5. Contact and Advanced Multibody — how Drake handles contact-rich systems

This section reviews hydroelastic contact, contact results, contact visualization, and more advanced contact dynamics.

- [`underactuated/book/00-tutorials/hydroelastic_contact_basics.ipynb`](underactuated/book/00-tutorials/hydroelastic_contact_basics.ipynb)  
  The main introduction to Drake's hydroelastic contact model, contact results, and `ContactVisualizer`.

- [`underactuated/book/00-tutorials/hydroelastic_contact_nonconvex_mesh.ipynb`](underactuated/book/00-tutorials/hydroelastic_contact_nonconvex_mesh.ipynb)  
  Shows how to configure hydroelastic contact for nonconvex mesh geometry.

---

#### 6. AutoDiff and Scalar Types

This section reviews `AutoDiffXd`, scalar conversion, gradients through systems, and gradients through multibody computations.

- [`underactuated/book/00-tutorials/autodiff_basics.ipynb`](underactuated/book/00-tutorials/autodiff_basics.ipynb)  
  The main introduction to AutoDiff in Drake, including gradient extraction and AutoDiff through systems.

- [`underactuated/book/00-tutorials/multibody_plant_autodiff_mass.ipynb`](underactuated/book/00-tutorials/multibody_plant_autodiff_mass.ipynb)  
  Demonstrates AutoDiff through `MultibodyPlant` computations with respect to mass parameters.

---

#### 7. MathematicalProgram Basics and Solver Workflow

This section reviews decision variables, costs, constraints, solvers, debugging, and repeated-solve workflows.

- [`underactuated/book/00-tutorials/mathematical_program.ipynb`](underactuated/book/00-tutorials/mathematical_program.ipynb)  
  The main introduction to `MathematicalProgram`, including variables, constraints, costs, solvers, and solution inspection.

- [`underactuated/book/00-tutorials/linear_program.ipynb`](underactuated/book/00-tutorials/linear_program.ipynb)  
  Builds and solves a small linear program as a template for LP modeling in Drake.

- [`underactuated/book/00-tutorials/quadratic_program.ipynb`](underactuated/book/00-tutorials/quadratic_program.ipynb)  
  Builds and solves a small quadratic program as a template for QP modeling in Drake.

- [`underactuated/book/00-tutorials/nonlinear_program.ipynb`](underactuated/book/00-tutorials/nonlinear_program.ipynb)  
  Introduces nonlinear programming, custom costs, custom constraints, and initial guesses.

- [`underactuated/book/00-tutorials/debug_mathematical_program.ipynb`](underactuated/book/00-tutorials/debug_mathematical_program.ipynb)  
  Shows how to inspect infeasible constraints, solver output, callbacks, and named bindings.

- [`underactuated/book/00-tutorials/solver_parameters.ipynb`](underactuated/book/00-tutorials/solver_parameters.ipynb)  
  Shows how to configure solver options and solver-specific parameters.

- [`underactuated/book/00-tutorials/updating_costs_and_constraints.ipynb`](underactuated/book/00-tutorials/updating_costs_and_constraints.ipynb)  
  Useful for repeated optimization workflows where costs or constraints are updated between solves.

- [`underactuated/book/00-tutorials/custom_gradients.ipynb`](underactuated/book/00-tutorials/custom_gradients.ipynb)  
  Shows how to provide custom gradients for optimization costs and constraints.

- [`underactuated/book/00-tutorials/sum_of_squares_optimization.ipynb`](underactuated/book/00-tutorials/sum_of_squares_optimization.ipynb)  
  Sets up a sum-of-squares program for certifying polynomial nonnegativity in Drake.

- [`underactuated/book/00-tutorials/licensed_solvers_deepnote.ipynb`](underactuated/book/00-tutorials/licensed_solvers_deepnote.ipynb)  
  Covers setup notes for licensed solvers such as Mosek and Gurobi.

</details>

---

<details>
<summary><strong>Humanoid / Legged Locomotion</strong></summary>

- [`slip.ipynb`](underactuated/book/04-simple_legs/slip.ipynb) — Builds the SLIP template as a `LeafSystem` whose `namedview` state and `MakeWitnessFunction` touchdown/takeoff/apex events encode the stance/flight hybrid switching behind the apex-to-apex return map.
- [`one_d_hopper.ipynb`](underactuated/book/04-simple_legs/exercises/one_d_hopper.ipynb) — Models the spring as a `LeafSystem` actuator and drives a `PreloadController` that detects bottom/apex from the body-velocity sign change and injects energy via a vectorized mechanical-energy budget.
- [`footstep_planning.ipynb`](underactuated/book/05-humanoids/exercises/footstep_planning.ipynb) — Decomposes MIQP footstep planning into stepping-stone halfspaces, one-hot stone binaries, big-M halfspace activation, step-span reachability limits, and a quadratic step-length cost solved by branch-and-bound.
- [`footstep_planning_gcs.ipynb`](underactuated/book/05-humanoids/exercises/footstep_planning_gcs.ipynb) — Recasts footstep planning as a GCS shortest path with `HPolyhedron` vertices, copied stone vertices for repeated steps, `e.xu()`/`e.xv()` edge reachability constraints, and unit edge costs under convex relaxation.
- [`littledog.ipynb`](underactuated/book/05-humanoids/littledog.ipynb) — A quadruped code study in generated `namedview` position/velocity views, per-gait `in_stance`/stride bookkeeping, per-timestep AutoDiff contexts, and whole-body (centroidal + full-kinematics) trajectory optimization with `PositionConstraint`/`OrientationConstraint`.
- [`compass_gait_limit_cycle.ipynb`](underactuated/book/17-contact/exercises/compass_gait_limit_cycle.ipynb) — Packages floating-base compass-gait dynamics into AutoDiff-compatible `MathematicalProgram` callbacks (manipulator equations, swing-foot kinematics, heel-strike impulse) with friction-cone contact forces and mirrored-periodicity constraints.
- [`basketball.ipynb`](underactuated/book/17-contact/basketball.ipynb) — A fixed-mode-sequence hybrid optimization that stitches analytic ballistic flight arcs together with `set_description`-labeled guard constraints and restitution/spin reset maps.
- [`multibody.ipynb`](underactuated/book/App-B-multibody/multibody.ipynb) — Compares time-stepping LCP contact resolution against MuJoCo-style relaxed complementarity-free contact by building the `q[n+1]` and contact-force `f[n]` update surfaces over the `(q, v)` grid.
- [`gcs.ipynb`](underactuated/book/App-C-optimization/gcs.ipynb) — A reusable Drake GCS pattern: `AddVertex` over `Point`/`VPolytope`/`Hyperellipsoid` sets, `AddEdge` costs, `e.xu()`/`e.xv()` edge variables, and a `SolveShortestPath` relaxation that reduces to the classic LP when every vertex is a point.

</details>

---

<details>
<summary><strong>Code Techniques Map</strong></summary>

Unique notebooks that are useful for code techs:

#### 1. Clean State Representation

How to make vector-valued robot states readable by using physically meaningful named fields instead of raw indices like `x[12]`.

- [`underactuated/book/04-simple_legs/slip.ipynb`](underactuated/book/04-simple_legs/slip.ipynb)  
  Uses `namedview` to access SLIP states as physical fields such as `s.r`, `s.theta`, and `s.rdot`.

- [`underactuated/book/10-trajopt/perching.ipynb`](underactuated/book/10-trajopt/perching.ipynb)  
  Uses named state views for glider dynamics and direct-collocation constraints, making trajectory constraints readable by physical coordinate.

- [`underactuated/book/05-humanoids/littledog.ipynb`](underactuated/book/05-humanoids/littledog.ipynb)  
  A rich example of generated named views for multibody positions, velocities, gait constraints, and stride construction.

---

#### 2. OOP and Problem Decomposition

How to separate physical data, model logic, optimization variables, and solver workflows into clean, readable components.

- [`underactuated/book/10-trajopt/exercises/orbital_transfer.ipynb`](underactuated/book/10-trajopt/exercises/orbital_transfer.ipynb)  
  Separates `Rocket`, `Planet`, and `Universe`, then uses clean dynamics and constraint residuals inside a trajectory optimization.

- [`underactuated/book/12-planning/exercises/rrt_planning.ipynb`](underactuated/book/12-planning/exercises/rrt_planning.ipynb)  
  Implements RRT/RRT* with clean algorithmic OOP, nested `Node` classes, and inheritance from `RRT` to `RRTStar`.

---

#### 3. Symbolic, Numeric, and AutoDiff Code Reuse

This section reviews how to write dynamics and constraint functions that can serve numerical rollout, symbolic derivative generation, and AutoDiff-based optimization.

- [`underactuated/book/10-trajopt/ilqr_cartpole.ipynb`](underactuated/book/10-trajopt/ilqr_cartpole.ipynb)  
  Writes cart-pole dynamics once and reuses it for numerical rollout and symbolic derivative generation with SymPy.

- [`underactuated/book/17-contact/exercises/compass_gait_limit_cycle.ipynb`](underactuated/book/17-contact/exercises/compass_gait_limit_cycle.ipynb)  
  Shows how to write multibody constraint callbacks that switch between `double` and `AutoDiffXd` plants for optimization.

- [`underactuated/book/10-trajopt/perching.ipynb`](underactuated/book/10-trajopt/perching.ipynb)  
  Uses a scalar-type-compatible custom glider plant, useful for trajectory optimization and feedback design workflows.

- [`underactuated/book/11-policy_search/policy_search.ipynb`](underactuated/book/11-policy_search/policy_search.ipynb)  
  Builds differentiable controller and running-cost systems, then differentiates through simulation rollouts using AutoDiff.

---

#### 4. Optimization Problem Architecture

This section reviews how to structure optimization code as a clear modeling pipeline: variables, constraints, costs, solve, and result extraction.

- [`underactuated/book/10-trajopt/exercises/orbital_transfer.ipynb`](underactuated/book/10-trajopt/exercises/orbital_transfer.ipynb)  
  Uses physical residual functions directly as trajectory-optimization constraints, keeping the dynamics and the constraint code in one place.

- [`underactuated/book/05-humanoids/exercises/footstep_planning.ipynb`](underactuated/book/05-humanoids/exercises/footstep_planning.ipynb)  
  Decomposes the MIQP footstep planner into a `SteppingStone`/`Terrain` geometry layer and named helper functions that add variables, reachability limits, big-M stone assignments, and costs.

- [`underactuated/book/03-acrobot/flatness.ipynb`](underactuated/book/03-acrobot/flatness.ipynb)  
  Builds a compact `PPTrajectory` wrapper that hides polynomial decision variables, continuity constraints, evaluation, and solving.

---

#### 5. Vectorization, Batching, and Neural Computation

This section reviews how to replace per-sample loops with batched arrays, `meshgrid`, matrix products, and `einsum`.

- [`underactuated/book/07-dp/exercises/pendulum_cvi.ipynb`](underactuated/book/07-dp/exercises/pendulum_cvi.ipynb)  
  Uses batched state grids, `BatchOutput`, and `np.einsum` to compute neural-network-style fitted value iteration efficiently.

- [`underactuated/book/07-dp/mlp.ipynb`](underactuated/book/07-dp/mlp.ipynb)  
  A companion fitted-value-iteration example using batched grids, vectorized costs, `BatchOutput`, backpropagation, and Adam.

- [`underactuated/book/06-stochastic/stochastic.ipynb`](underactuated/book/06-stochastic/stochastic.ipynb)  
  Simulates many stochastic particles as one vectorized Drake system, with packed particle states and a matching custom visualizer.

- [`underactuated/book/08-lqr/value_iteration.ipynb`](underactuated/book/08-lqr/value_iteration.ipynb)  
  A compact vectorized value-iteration reference using `meshgrid` and vectorized quadratic forms.

</details>

---

<details>
<summary><strong>Some useful links</strong></summary>

https://underactuated.csail.mit.edu/

https://github.com/RussTedrake/underactuated

https://underactuated.csail.mit.edu/Spring2024/

https://www.youtube.com/channel/UChfUOAhz7ynELF-s_1LPpWg

</details>
