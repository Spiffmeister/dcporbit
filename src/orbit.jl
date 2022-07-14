module orbit

    using LinearAlgebra
    using Roots
    using Plots
    using LaTeXStrings
    using JLD2
    using Distributed

    include("Particle.jl")
    include("analytic.jl")
    include("event_location.jl")
    include("Integrators.jl")
    include("Plotting.jl")

    export sim,particle,exact_particle,forces
    export analytic_solve,integrate!,solve_orbit!,run_sim!
    export MagneticForce,MagneticForce_GC

end