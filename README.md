# SEMDynamics.jl

| Documentation | Build | License |
|:---:|:---:|:---:|
| [![Stable documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://vanyansin2000.github.io/SEMDynamics.jl/stable/) [![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://vanyansin2000.github.io/SEMDynamics.jl/dev/) | [![CI](https://github.com/vanyansin2000/SEMDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vanyansin2000/SEMDynamics.jl/actions/workflows/CI.yml) | [![License: MPL-2.0](https://img.shields.io/badge/license-MPL--2.0-green.svg)](LICENSE) |

`SEMDynamics.jl` is a Julia package for numerical trajectory analysis in the
Earth--Moon restricted multibody problem.  It provides equations of motion for
the circular restricted three-body problem (CR3BP) and the Earth--Moon
bicircular restricted four-body problem (BCR4BP), together with state
transition-matrix propagation, event callbacks, coordinate transformations,
trajectory plotting, and differential correction of planar distant retrograde
orbits (DROs).

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
uf , ts , us = integration(u0, (0.0, 10.0), ode_params(cr3bp_eqm! , ode_args , Bcr4bp_aux())) 
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
more corrections.  A closer two-element seed `[x, vy]` and its corresponding
period can also be supplied explicitly:

```julia
orbit = generate_DRO(
    P=3.4,
    seed=[1.1754, -0.4943],
    seed_period=pi,
)
```

The current high-level generator is limited to planar CR3BP DROs, whose initial
state is `[x, 0, 0, 0, vy, 0]`.  The lower-level
`orbit_shooting(residual!, initial_guess, parameters)` interface can be used to
define differential corrections for other orbit families.

### Detect events

The package supplies callbacks compatible with DifferentialEquations.jl.

```julia
using SEMDynamics
using DifferentialEquations

aux = Bcr4bp_Aux()
events = dynamic_events()
callback = CallbackSet(cb_p2collision(events), cb_escape(events))

u0 = [1.05, 0.0, 0.0, 0.2]
problem = ODEProblem(bcr4bp_eqm!, u0, (0.0, 10.0), aux)
solution = solve(problem, Vern7(); callback, reltol=1e-12, abstol=1e-12)
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
| `orbit_shooting(...)` | General nonlinear shooting interface. |
| `compute_jacobi` | Computes the CR3BP Jacobi constant. |
| `cb_*` | Collision, escape, energy-crossing, and apsis event callbacks. |
| `cr3bp_inertial_to_rotating` / `cr3bp_rotating_to_inertial` | Converts CR3BP states between inertial and rotating frames. |
| `plotPO!`, `plot_traj_rot!`, `plot_traj_inertial!` | Makie plotting helpers for periodic or integrated trajectories. |

Use Julia help mode to inspect exported methods and their arguments:

```julia
help?> generate_DRO
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
