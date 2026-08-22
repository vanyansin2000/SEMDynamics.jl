"""
    SEMDynamics

Earth--Moon CR3BP/BCR4BP dynamics, event detection, coordinate transforms,
periodic-orbit generation, numerical utilities, and Makie visualization.
"""
module SEMDynamics

using Reexport

include("Utils.jl")
include("Dynamics.jl")
include("PeriodicOrbits.jl")
include("Visualization.jl")


@reexport using .Utils
@reexport using .Dynamics
@reexport using .PeriodicOrbits
@reexport using .Visualization

end
