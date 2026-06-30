## Underactuated Robotics

This repository collects my worked notebooks, implementation notes, and small companion code for studying Underactuated Robotics. Material is organized under [`underactuated/book/`](underactuated/book/).

---

### Map

Repo notebook map for [`underactuated/book/`](underactuated/book/).

```text
underactuated/book/
|-- 00-tutorials/ - Drake tutorials; https://underactuated.csail.mit.edu/drake.html
|  |-- index.ipynb - tutorial index; Drake tutorial index.
|  |-- dynamical_systems.ipynb - Drake tutorial; LeafSystem and simulation basics.
|  |-- authoring_leaf_systems.ipynb - Drake tutorial; Custom LeafSystem interface.
|  |-- working_with_diagrams.ipynb - Drake tutorial; Diagram and Context composition.
|  |-- authoring_multibody_simulation.ipynb - Drake tutorial; MultibodyPlant simulation template.
|  |-- pyplot_animation_multibody_plant.ipynb - Drake tutorial; PyPlot multibody animation.
|  |-- rendering_multibody_plant.ipynb - Drake tutorial; RGB-D and label rendering.
|  |-- configuring_rendering_lighting.ipynb - Drake tutorial; PBR lighting configuration.
|  |-- hydroelastic_contact_basics.ipynb - Drake tutorial; Hydroelastic contact basics.
|  |-- hydroelastic_contact_nonconvex_mesh.ipynb - Drake tutorial; Nonconvex mesh contact.
|  |-- autodiff_basics.ipynb - Drake tutorial; AutoDiff basics.
|  |-- multibody_plant_autodiff_mass.ipynb - Drake tutorial; Mass-parameter AutoDiff.
|  |-- mathematical_program.ipynb - Drake tutorial; MathematicalProgram core.
|  |-- linear_program.ipynb - Drake tutorial; LP modeling template.
|  |-- quadratic_program.ipynb - Drake tutorial; QP modeling template.
|  |-- nonlinear_program.ipynb - Drake tutorial; NLP custom constraints.
|  |-- sum_of_squares_optimization.ipynb - Drake tutorial; SOS/SDP basics.
|  |-- debug_mathematical_program.ipynb - Drake tutorial; Optimization debugging.
|  |-- solver_parameters.ipynb - Drake tutorial; Solver option setup.
|  |-- updating_costs_and_constraints.ipynb - Drake tutorial; Update bindings between solves.
|  |-- custom_gradients.ipynb - Drake tutorial; Manual optimization gradients.
|  |-- mathematical_program_multibody_plant.ipynb - Drake tutorial; IIWA IK with AutoDiff evaluator.
|  `-- licensed_solvers_deepnote.ipynb - Drake tutorial; Commercial solvers on Deepnote.

|-- 01-intro/ - Chapter 1, Fully-actuated vs Underactuated Systems; https://underactuated.csail.mit.edu/intro.html
|  |-- intro.ipynb - notes example; Double pendulum dynamics intro.
|  `-- exercises/drake_systems.ipynb - Exercise 1.5 / Spring set 1; Write a simple Drake System.

|-- 02-pend/ - Chapter 2, The Simple Pendulum; https://underactuated.csail.mit.edu/pend.html
|  |-- pend.ipynb - notes example; Simple pendulum simulation.
|  |-- attractivity.ipynb - notes example; Attractive unstable equilibrium.
|  |-- autapse.ipynb - notes example; Single-neuron autapse dynamics.
|  |-- lstm.ipynb - notes example; LSTM/JANET dynamics view.
|  |-- energy_shaping.ipynb - notes example; Pendulum energy shaping swing-up.
|  |-- exercises/attractivity_vs_stability.ipynb - Exercise 2.4; Attractivity vs stability.
|  |-- exercises/vibrating_pendulum.ipynb - Exercise 2.5 / Spring set 1; Vibrating-base pendulum control.
|  `-- exercises/hopfield_network.ipynb - Exercise 2.6; Hopfield network image recovery.

|-- 03-acrobot/ - Chapter 3, Acrobots, Cart-Poles, and Quadrotors; https://underactuated.csail.mit.edu/acrobot.html
|  |-- acrobot.ipynb - notes example; Acrobot dynamics + LQR.
|  |-- cartpole.ipynb - notes example; Cart-Pole dynamics + Meshcat.
|  |-- cartpole_energy_shaping.ipynb - local example; Cart-Pole energy shaping.
|  |-- planar_quadrotor.ipynb - notes example; Planar quadrotor LQR.
|  |-- quadrotor.ipynb - notes example; 3D quadrotor MultibodyPlant.
|  |-- flatness.ipynb - notes example; Differential flatness template.
|  |-- exercises/cartpole_balancing.ipynb - Exercise 3.1 / Spring set 3; Cart-Pole LQR balancing.
|  `-- exercises/cartpoles_urdf.ipynb - Exercise 3.2 / Spring set 3; Write a cart-pole URDF.

|-- 04-simple_legs/ - Chapter 4, Simple Models of Walking and Running; https://underactuated.csail.mit.edu/simple_legs.html
|  |-- rimless_wheel.ipynb - notes example; Rimless wheel limit cycle.
|  |-- compass_gait.ipynb - notes example; Passive compass gait simulation.
|  |-- slip.ipynb - notes example; SLIP hybrid events + map.
|  |-- planar_one_leg_hopper.ipynb - local example; Planar hopper control draft.
|  `-- exercises/one_d_hopper.ipynb - Exercise 4.1 / Spring set 8; Raibert 1D hopper control.

|-- 05-humanoids/ - Chapter 5, Highly-articulated Legged Robots; https://underactuated.csail.mit.edu/humanoids.html
|  |-- zmp_planner.ipynb - notes example; ZMP CoM trajectory planning.
|  |-- littledog.ipynb - notes example; LittleDog gait optimization.
|  |-- spot.ipynb - local example; Spot model visualization entry.
|  |-- exercises/footstep_planning.ipynb - Exercise 5.1; MIQP footstep planning.
|  `-- exercises/footstep_planning_gcs.ipynb - Exercise 5.3 / Spring set 10; GCS footstep planning.

|-- 06-stochastic/ - Chapter 6, Model Systems with Stochasticity; https://underactuated.csail.mit.edu/stochastic.html
|  `-- stochastic.ipynb - notes example; Batch stochastic particle simulation.

|-- 07-dp/ - Chapter 7, Dynamic Programming; https://underactuated.csail.mit.edu/dp.html
|  |-- grid_world.ipynb - notes example; Grid-world value iteration.
|  |-- on_a_mesh.ipynb - notes example; Mesh-based continuous DP.
|  |-- mlp.ipynb - notes example; Neural fitted value iteration.
|  |-- exercises/minimum_time.ipynb - Exercise 7.5 / Spring set 2; Minimum-time value iteration.
|  |-- exercises/lp_dp.ipynb - Exercise 7.6 / Spring set 2; DP as a linear program.
|  `-- exercises/pendulum_cvi.ipynb - Exercise 7.7 / Spring set 4; Pendulum continuous FVI.

|-- 08-lqr/ - Chapter 8, Linear Quadratic Regulators; https://underactuated.csail.mit.edu/lqr.html
|  |-- continuous_vs_discrete_time.ipynb - notes example; Continuous/discrete LQR comparison.
|  |-- value_iteration.ipynb - notes example; FVI reproduces LQR.
|  |-- ballbot.ipynb - notes example; Ballbot reduced LQR.
|  |-- manifold.ipynb - local example; LQR on manifold coordinates.
|  `-- exercises/drake_diagrams.ipynb - Exercise 8.1 / Spring set 2; Drake Diagram + LQR.

|-- 09-lyapunov/ - Chapter 9, Lyapunov Analysis; https://underactuated.csail.mit.edu/lyapunov.html
|  |-- common_lyap_linear.ipynb - notes example; Common quadratic Lyapunov.
|  |-- global_polynomial.ipynb - notes example; SOS Lyapunov search.
|  |-- cubic_poly.ipynb - notes example; Cubic ROA certificate.
|  |-- van_der_pol_w_alternations.ipynb - notes example; Van der Pol alternating SOS.
|  |-- star_convex.ipynb - notes example; Nonconvex ROA example.
|  |-- outer_approx.ipynb - notes example; Convex outer approximation.
|  |-- global_pend.ipynb - notes example; Pendulum global stability.
|  |-- approximate_dp.ipynb - notes example; SOS approximate DP.
|  |-- sampling.ipynb - local example; Sampling + LP stability.
|  |-- trig_poly.ipynb - local example; Trig polynomial verification.
|  |-- exercises/control.ipynb - Exercise 9.8; Wheeled robot CLF control.
|  |-- exercises/sos_and_psd.ipynb - Exercise 9.9 / Spring set 6; PSD vs SOS polynomials.
|  `-- exercises/van_der_pol.ipynb - Exercise 9.10 / Spring set 5; Van der Pol SOS ROA.

|-- 10-trajopt/ - Chapter 10, Trajectory Optimization; https://underactuated.csail.mit.edu/trajopt.html
|  |-- double_integrator.ipynb - notes example; Direct transcription/collocation intro.
|  |-- dircol.ipynb - notes example; Drake dircol swing-up template.
|  |-- perching.ipynb - notes example; Glider trajopt + finite-horizon LQR.
|  |-- perching/perching.ipynb - notes example; Perching companion notebook.
|  |-- mi_convex.ipynb - notes example; Mixed-integer collision avoidance.
|  |-- gcs_quadrotor.ipynb - notes example; Quadrotor GCS planning + rounding.
|  |-- ilqr_cartpole.ipynb - local example; Cart-Pole iLQR swing-up.
|  |-- cartpole/MATLAB/ilqr_cartpole.ipynb - local example; MATLAB-style iLQR Cart-Pole.
|  |-- cartpole/ct-scvx.ipynb - local example; Continuous-time SCvx.
|  |-- compare_dirtrans_dircol.ipynb - local example; Compare transcription and collocation.
|  |-- exercises/shooting_vs_transcription.ipynb - Exercise 10.1; Compare shooting/transcription.
|  |-- exercises/orbital_transfer.ipynb - Exercise 10.2 / Spring set 6; Earth-to-Mars rocket trajopt.
|  `-- exercises/ilqr_driving.ipynb - Exercise 10.4 / Spring set 7; Autonomous driving iLQR.

|-- 11-policy_search/ - Chapter 11, Policy Search; https://underactuated.csail.mit.edu/policy_search.html
|  `-- policy_search.ipynb - notes example; LQR policy optimization.

|-- 12-planning/ - Chapter 12, Sampling-based motion planning; https://underactuated.csail.mit.edu/planning.html
|  `-- exercises/rrt_planning.ipynb - Exercise 12.1; Implement RRT/RRT*.

|-- 13-robust/ - Chapter 13, Robust and Stochastic Control; https://underactuated.csail.mit.edu/robust.html
|  `-- quadrotor_in_wind.ipynb - notes example; Wind disturbance LQR.

|-- 15-output_feedback/ - Chapter 15, Output Feedback (aka Pixels-to-Torques); https://underactuated.csail.mit.edu/output_feedback.html
|  `-- acrobot_w_encoders.ipynb - local example; Encoder-only Acrobot.

|-- 16-limit_cycles/ - Chapter 16, Algorithms for Limit Cycles; https://underactuated.csail.mit.edu/limit_cycles.html
|  `-- limit_cycles.ipynb - notes example; Van der Pol limit cycle trajopt.

|-- 17-contact/ - Chapter 17, Planning and Control through Contact; https://underactuated.csail.mit.edu/contact.html
|  |-- rimless_wheel.ipynb - notes example; Contact-aware rimless trajopt.
|  |-- basketball.ipynb - notes example; Basketball hybrid trajopt.
|  |-- hybrid.ipynb - local example; Hybrid multibody collocation.
|  `-- exercises/compass_gait_limit_cycle.ipynb - Exercise 17.1 / Spring set 8; Compass gait limit-cycle NLP.

|-- 18-sysid/ - Chapter 18, System Identification; https://underactuated.csail.mit.edu/sysid.html
|  |-- sysid.ipynb - notes example; Acrobot/Cart-Pole sysid.
|  |-- exercises/linear_sysid.ipynb - Exercise 18.1; Linear A/B least-squares.
|  `-- exercises/glider_sysid.ipynb - Exercise 18.2; Glider basis-function sysid.

|-- App-B-multibody/ - Appendix B, Multi-Body Dynamics; https://underactuated.csail.mit.edu/multibody.html
|  `-- multibody.ipynb - notes example; LCP vs relaxed contact.

|-- App-C-optimization/ - Appendix C, Optimization and Mathematical Programming; https://underactuated.csail.mit.edu/optimization.html
|  |-- sdp.ipynb - notes example; SDP relaxation template.
|  |-- sos_six_hump_camel.ipynb - notes example; Six-Hump Camel via SOS.
|  `-- gcs.ipynb - notes example; GCS shortest path template.

`-- figures/ - figure helper notebooks; https://underactuated.csail.mit.edu/
   |-- Quadrotor2D.ipynb - figure helper; Planar quadrotor figure.
   `-- lcp_cart.ipynb - figure helper; LCP cart contact surface.```

---

<details>
<summary><strong>Links</strong></summary>

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
<summary><strong>Humanoid</strong></summary>

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
<summary><strong>Code Tech Map</strong></summary>

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
<summary><strong>Hydroelastic Contact</strong></summary>

Hydroelastic contact models contact as a pressure-distributed surface patch rather than a small set of point contacts.

- [Newton documentation: Collisions and Contacts](https://newton-physics.github.io/newton/stable/concepts/collisions.html)  
  The best Newton-specific starting point. It explains how hydroelastic contact is enabled in Newton, including SDF requirements, mesh/primitive setup, contact flags, stiffness parameters, and example scripts.

- [GPU-Accelerated Hydroelastic Contact via Signed Distance Fields](https://openreview.net/pdf?id=ogndqznZyY)  
  The key Newton-related technical paper. It explains Newton's SDF-based hydroelastic contact formulation, including how contact surfaces are extracted and reduced for GPU-parallel simulation.

- [NVIDIA Technical Blog: Newton Adds Contact-Rich Manipulation and Locomotion Capabilities](https://developer.nvidia.com/blog/newton-adds-contact-rich-manipulation-and-locomotion-capabilities-for-industrial-robotics/)  
  A practical engineering overview of Newton 1.0, SDF collision, hydroelastic contact, MuJoCo Warp integration, and Isaac Lab usage.

- [TRI Medium: Rethinking Contact Simulation for Robot Manipulation](https://medium.com/toyotaresearch/rethinking-contact-simulation-for-robot-manipulation-434a56b5ec88)  
  The most intuitive introduction to why point contact can be insufficient for manipulation and why contact patches and pressure fields are useful.

- [Drake Hydroelastic Contact User Guide](https://drake.mit.edu/doxygen_cxx/group__hydroelastic__user__guide.html)  
  A strong implementation-oriented reference for the original Drake-style hydroelastic contact model, including compliant/rigid hydroelastic geometry choices and fallback behavior.

- [Drake hydroelastic contact tutorial notebook](https://github.com/RobotLocomotion/drake/blob/master/tutorials/hydroelastic_contact_basics.ipynb)  
  A hands-on tutorial for setting up simulations with hydroelastic contact and inspecting contact results.

- [A Pressure Field Model for Fast, Robust Approximation of Net Contact Force and Moment Between Nominally Rigid Objects](https://arxiv.org/abs/1904.11433)  
  The foundational hydroelastic / pressure-field contact paper. It introduces the idea of object-centric pressure fields and contact surfaces defined by equal pressure.

- [Velocity Level Approximation of Pressure Field Contact Patches](https://arxiv.org/abs/2110.04157)  
  An important follow-up that makes pressure-field contact compatible with velocity-level time-stepping solvers, which is crucial for practical multibody simulation.

</details>

---

<details>
<summary><strong>Official links</strong></summary>

https://underactuated.csail.mit.edu/

https://github.com/RussTedrake/underactuated

https://underactuated.csail.mit.edu/Spring2024/

https://openreview.net/pdf?id=ogndqznZyY

https://www.youtube.com/channel/UChfUOAhz7ynELF-s_1LPpWg

</details>

</details>

---