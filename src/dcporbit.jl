"""
    dcporbit

Module for solving the particle orbit
"""
module dcporbit

    using LinearAlgebra
    using Roots
    using Plots
    using LaTeXStrings
    using JLD2
    using Distributed

    include("Particle.jl")
    include("Analytic.jl")
    include("event_location.jl")
    include("Integrators.jl")
    include("timesolver.jl")
    include("Plotting.jl")

    export FullOrbit, FO, GuidingCentre, GC
    export sim,particle,exact_particle,forces
    export analytic_solve,integrate!,solve_orbit!,run_sim!
    export MagneticForce,MagneticForce_GC

end