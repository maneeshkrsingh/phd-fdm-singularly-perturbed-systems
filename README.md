# PhD FDM Code: Singularly Perturbed Systems

MATLAB code developed during PhD research on **Finite Difference Methods (FDM)** for
singularly perturbed systems of convection-diffusion and reaction-diffusion equations.

## Repository Structure

### `convection-diffusion/`
FDM schemes for singularly perturbed convection-diffusion systems.

| Subfolder | Description |
|---|---|
| `basic-discretizations/` | Upwind and basic schemes on uniform/piecewise-uniform meshes (Q1–Q4 programs, discontinuous source) |
| `shishkin-mesh/` | Upwind FDM on Shishkin and Bakhvalov-Shishkin (BS) meshes; single/double/triple perturbation parameters |
| `hybrid/` | Hybrid differencing schemes (midpoint upwind + central difference) on Shishkin and BS meshes |
| `richardson-extrapolation/` | Richardson extrapolation to improve order of convergence for convection-diffusion systems |
| `discontinuous-source/` | FDM for systems with discontinuous convection coefficients / source terms (time-dependent) |
| `2d-fractional-step/` | Fractional-step (ADI-type) methods for 2D singularly perturbed convection-diffusion systems |
| `semilinear/` | Semilinear singularly perturbed convection-diffusion systems on Shishkin and BS meshes |
| `ode-trial/` | ODE solver-based (MATLAB `pde/ode` toolbox) trials for convection-diffusion systems |
| `equidistribution/` | Equidistribution mesh generation for convection-diffusion problems |
| `delay-richardson/` | Richardson extrapolation for delay convection-diffusion equations |

### `reaction-diffusion/`
FDM schemes for singularly perturbed reaction-diffusion systems.

| Subfolder | Description |
|---|---|
| `basic/` | Basic upwind/central FDM for 1D reaction-diffusion systems (smooth and non-smooth cases) |
| `2d-rd-fully-discrete/` | Fully discrete (space-time) FDM for 2D reaction-diffusion systems |
| `2d-elliptic/` | 2D elliptic singularly perturbed reaction-diffusion systems; fractional-step methods |
| `2d-fractional-step/` | Fractional-step ADI-type methods for 2D reaction-diffusion systems |
| `gracia-lisbona-system/` | Numerical experiments for Gracia-Lisbona coupled reaction-diffusion systems |

### `discontinuous-galerkin/`
Finite Element Method (FEM) and Discontinuous Galerkin (DG) trials.

| Subfolder | Description |
|---|---|
| `dg-trial/` | Simple DG solver trials |
| `gautam-dgfem/` | DG-FEM implementation for coupled convection-reaction systems |
| `fem-system/` | 1D FEM for coupled singularly perturbed systems (coupled, scalar, reaction subproblems) |
| `fem-trial/` | FEM trial programs for convection-reaction SPPs |
| `trial-fem/` | Additional FEM trials (convection, reaction, general FEM code) |

### `basic-schemes/`
Foundational FDM schemes and utilities.

| Subfolder / File | Description |
|---|---|
| `*.m` (root) | FTCS, FTBS, BTCS, Crank-Nicolson, upwind-fitted, Lax-Wendroff, Shishkin degenerate, etc. |
| `rao-hodie/` | Rao HODIE (High Order Difference with Identity Expansion) and Richardson-accelerated versions |
| `equidistribution/` | Equidistribution-based adaptive mesh methods |
| `equidistribution-variants/` | Alternative equidistribution implementations |
| `richardson/` | Richardson extrapolation on Shishkin mesh for scalar SPPs |
| `experiments/` | Experimental / comparative studies of various schemes |

### `utils/`
Utility and helper scripts.

| File/Folder | Description |
|---|---|
| `sys-ode/` | ODE system helper functions for coupled convection-diffusion |
| `jusr.m` | Miscellaneous utility |

## Key Topics Covered

- Singularly perturbed systems with **single, double, and triple perturbation parameters**
- **Shishkin, Bakhvalov-Shishkin (BS), modified BS (MBS), and Vulanovic-Shishkin (VS)** meshes
- **Upwind, hybrid, and HODIE** finite difference operators
- **Richardson extrapolation** for higher-order accuracy
- **Discontinuous convection coefficients and source terms** (interior layers)
- **Time-dependent** (parabolic) and **steady-state** (elliptic) formulations
- **2D problems** via fractional-step (ADI) and fully-discrete methods
- **FEM and DG** methods as complementary approaches

## Requirements

- MATLAB R2014b or later

## Usage

Each subfolder is self-contained. Open the main `.m` file in MATLAB and run directly.
Mesh generation functions (`BSmesh.m`, `VSmesh.m`, `MBSmesh.m`, etc.) are co-located
with the solver programs that use them.
