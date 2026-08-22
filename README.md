# SEMDynamics.jl

English | [中文](README.zh.md)

| Documentation | Build | License |
|:---:|:---:|:---:|
| [![Stable documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://vanyansin2000.github.io/SEMDynamics.jl/stable/) [![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://vanyansin2000.github.io/SEMDynamics.jl/dev/) | [![CI](https://github.com/vanyansin2000/SEMDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vanyansin2000/SEMDynamics.jl/actions/workflows/CI.yml) | [![License: MPL-2.0](https://img.shields.io/badge/license-MPL--2.0-green.svg)](LICENSE) |

`SEMDynamics.jl` is a Julia package for numerical trajectory analysis in the
Earth--Moon restricted multibody problem.  It provides equations of motion for
the circular restricted three-body problem (CR3BP) and the Earth--Moon
bicircular restricted four-body problem (BCR4BP), together with state
transition-matrix propagation, event callbacks, coordinate transformations,
trajectory plotting, and differential correction of planar distant retrograde
orbits (DROs) and spatial halo orbits.

Unless stated otherwise, states are expressed in the nondimensional rotating
frame and times and periods are nondimensional.

## Installation

After the package has been registered, install it from the Julia General
registry:

```julia
using Pkg
Pkg.add("SEMDynamics")
```

SEMDynamics requires Julia 1.12 or later.

## Quick start

Integrate a planar CR3BP trajectory 

```julia
using SEMDynamics

u0 = [0.8, 0.1, 0.02, -0.03] # [x, y, vx, vy]
ode_args = (; reltol = 1e-12 , abstol = 1e-12) # default 
uf , ts , us = integration(u0, (0.0, 10.0), ode_params(cr3bp_eqm! , ode_args , Bcr4bp_Aux()))
# or simply
uf , ts , us = integration(u0, (0.0, 10.0), ode_params(cr3bp_eqm!)) 
```

### Generate a planar DRO

`generate_DRO` differentially corrects a planar CR3BP distant retrograde orbit
with the requested full period.  The default seed is a 2:1 DRO with period
`pi`.

```julia
using SEMDynamics

orbit = generate_DRO(P=pi)
orbit.x0 # six-component initial state in the rotating frame
orbit.P  # nondimensional full period
orbit.C  # Jacobi constant
```

For a target period away from the reference seed, the solver performs
short-step continuation in period, using a secant predictor followed by Newton
correction at every step:

```julia
orbit = generate_DRO(P=3.4, continuation_step=0.05)
```

Smaller values of `continuation_step` are generally more robust but require
more corrections. Every periodic-orbit seed uses the same complete six-state
form `[x, y, z, vx, vy, vz]`; a closer seed and its corresponding period can
also be supplied explicitly:

```julia
orbit = generate_DRO(
    P=3.4,
    seed=[1.1754, 0.0, 0.0, 0.0, -0.4943, 0.0],
    seed_period=pi,
)
```

### Generate halo orbits and a 9:2 NRHO

Halo families are selected by their `branch` and libration point. The initial
state is a crossing of the x-z symmetry plane, `[x, 0, z, 0, vy, 0]`.

```julia
northern_l1 = generate_halo(branch=:northern, lp=:L1)
southern_l2 = generate_halo(branch=:southern, lp=:L2)

# Continue the selected family to another nondimensional period.
continued = generate_halo(branch=:northern, lp=:L1, P=2.75)

# Continue the southern L2 family to its 9:2 lunar-synodic resonant period.
nrho = generate_nrho_9_2()
```

Halo and DRO differential correction use the same unknowns `[x, z, vy]`, the
same half-period conditions `[y, vx, vz] = 0`, and a default shooting tolerance
of `1e-10`. Halo period continuation defaults to a step of `0.005`.

The lower-level `orbit_shooting(dynamics, state_guess, P)` interface accepts a
complete six-state and returns the corrected symmetry-plane state. The generic
`orbit_shooting(residual!, initial_guess, parameters)` overload remains
available for custom corrections.

### Detect events

The package supplies callbacks compatible with DifferentialEquations.jl.

```julia
using SEMDynamics
using DifferentialEquations

aux = Bcr4bp_Aux()
events = dynamic_events()
callback = CallbackSet(
    cb_p2collision(events),
    cb_enter(events; terminate=false),
    cb_escape(events),
)

u0 = [1.05, 0.0, 0.0, 0.2]
problem = ODEProblem(bcr4bp_eqm!, u0, (0.0, 10.0), aux)
solution = solve(problem, Vern7(); callback, reltol=1e-12, abstol=1e-12)
```

`cb_enter` and `cb_escape` detect inward and outward crossings of the same P2
sphere. Lunar apsides use the same combined interface as terrestrial apsides:

```julia
earth_apses = cb_apse_p1(events; terminate_perigee=false, terminate_apogee=false)
lunar_apses = cb_apse_p2(events; terminate_perilune=false, terminate_apolune=false)
```

### Plot an orbit

```julia
using CairoMakie
using SEMDynamics

orbit = generate_DRO(P=pi)
figure = Figure()
axis = Axis(figure[1, 1], aspect=AxisAspect(1))
cr3bp_set_scn!(axis)
plotPO!(axis, orbit; color=:crimson, linewidth=2)
display(figure)
```

## State conventions

| Model | State | State with STM |
| --- | --- | --- |
| Planar CR3BP/BCR4BP | `[x, y, vx, vy]` | 20 elements: state followed by a 4x4 STM |
| Spatial CR3BP/BCR4BP | `[x, y, z, vx, vy, vz]` | 42 elements: state followed by a 6x6 STM |

`cr3bp_eqm!` and `bcr4bp_eqm!` dispatch automatically from the state length.
For the CR3BP, `compute_jacobi` evaluates the Jacobi constant.

## Main interface

| Interface | Description |
| --- | --- |
| `Bcr4bp_Aux()` | Creates physical constants and nondimensional parameters for the Earth--Moon model. |
| `cr3bp_eqm!` | CR3BP equations of motion for planar, spatial, and STM-augmented states. |
| `bcr4bp_eqm!` | Earth--Moon BCR4BP equations of motion for planar, spatial, and STM-augmented states. |
| `generate_DRO(P=...)` | Generates a differentially corrected planar CR3BP DRO. |
| `generate_halo(branch=..., lp=...)` | Generates northern or southern L1/L2 CR3BP halo orbits. |
| `generate_nrho_9_2()` | Continues the L2 halo branch to its 9:2 lunar-synodic period. |
| `orbit_shooting(...)` | General nonlinear shooting interface. |
| `compute_jacobi` | Computes the CR3BP Jacobi constant. |
| `cb_enter` / `cb_escape` | Detects inward/outward crossings of the P2 sphere. |
| `cb_apse_p1` / `cb_apse_p2` | Detects perigee/apogee or perilune/apolune events. |
| `cb_*` | Collision, energy-crossing, section-crossing, and other event callbacks. |
| `cr3bp_inertial_to_rotating` / `cr3bp_rotating_to_inertial` | Converts CR3BP states between inertial and rotating frames. |
| `plotPO!`, `plot_traj_rot!`, `plot_traj_inertial!` | Makie plotting helpers for periodic or integrated trajectories. |

Use Julia help mode to inspect exported methods and their arguments:

```julia
help?> generate_DRO
help?> generate_halo
```

## Documentation

The online documentation is available at
[vanyansin2000.github.io/SEMDynamics.jl](https://vanyansin2000.github.io/SEMDynamics.jl/stable/).

To build the documentation locally from the repository root:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated site is written to `docs/build/`.

## Testing

Run the package test suite from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Contributing

Issues and pull requests are welcome.  Please include tests for behavior
changes and document public APIs with Julia docstrings.  
<!-- See [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) for the maintainer release and
registration checklist. -->

## License

SEMDynamics.jl is distributed under the [Mozilla Public License 2.0](LICENSE).
