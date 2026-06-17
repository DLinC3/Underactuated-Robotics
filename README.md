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

- [`underactuated/book/App-B-multibody/`](underactuated/book/App-B-multibody/)  
  A deeper multibody appendix for dynamics, contact, and complementarity-related ideas.

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

- [`underactuated/book/17-contact/`](underactuated/book/17-contact/)  
  A broader contact chapter covering hybrid dynamics, contact examples, and optimization-based contact problems.

- [`underactuated/book/App-B-multibody/`](underactuated/book/App-B-multibody/)  
  A deeper appendix for multibody dynamics and contact-related theory.

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
  Introduces linear programming in Drake.

- [`underactuated/book/00-tutorials/quadratic_program.ipynb`](underactuated/book/00-tutorials/quadratic_program.ipynb)  
  Introduces quadratic programming in Drake.

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
  Introduces sum-of-squares optimization in Drake.

- [`underactuated/book/00-tutorials/licensed_solvers_deepnote.ipynb`](underactuated/book/00-tutorials/licensed_solvers_deepnote.ipynb)  
  Covers setup notes for licensed solvers such as Mosek and Gurobi.

- [`underactuated/book/App-c-optimization/`](underactuated/book/App-c-optimization/)  
  A deeper optimization appendix covering advanced topics beyond the basic `MathematicalProgram` tutorials.


</details>

---

<details>
<summary><strong>Code Tech Map</strong></summary>

Unique notebooks that are useful for code techs:

#### 1. Clean State Representation

How to make vector-valued robot states readable by using physically meaningful named fields instead of like x[12].

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
  Implements RRT/RRT* with clean algorithmic OOP, nested `Node` classes, and inheritance from `RRT` to `RRTStar` (I feel classic but useful). 

- [`underactuated/book/05-humanoids/exercises/footstep_planning.ipynb`](underactuated/book/05-humanoids/exercises/footstep_planning.ipynb)  
  Separates terrain geometry, stepping-stone data, and mixed-integer optimization constraints into readable layers.

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

#### 3. Optimization Problem Architecture

This section reviews how to structure optimization code as a clear modeling pipeline: variables, constraints, costs, solve, and result extraction.

- [`underactuated/book/10-trajopt/exercises/orbital_transfer.ipynb`](underactuated/book/10-trajopt/exercises/orbital_transfer.ipynb)  
  A clean example of using physical residual functions directly as trajectory-optimization constraints.

- [`underactuated/book/05-humanoids/exercises/footstep_planning.ipynb`](underactuated/book/05-humanoids/exercises/footstep_planning.ipynb)  
  Shows how to decompose a mixed-integer footstep planner into geometry classes and named optimization helper functions.

- [`underactuated/book/17-contact/exercises/compass_gait_limit_cycle.ipynb`](underactuated/book/17-contact/exercises/compass_gait_limit_cycle.ipynb)  
  Demonstrates careful vector packing, AutoDiff-compatible constraints, and multibody callbacks inside `MathematicalProgram`.

- [`underactuated/book/03-acrobot/flatness.ipynb`](underactuated/book/03-acrobot/flatness.ipynb)  
  Defines a compact `PPTrajectory` wrapper that hides polynomial decision variables, continuity constraints, evaluation, and solving.

---

#### 4. Vectorization, Batching, and Neural Computation

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
<summary><strong>underactuated links</strong></summary>

https://underactuated.csail.mit.edu/

https://github.com/RussTedrake/underactuated

https://underactuated.csail.mit.edu/Spring2024/

https://openreview.net/pdf?id=ogndqznZyY

https://www.youtube.com/channel/UChfUOAhz7ynELF-s_1LPpWg

</details>
