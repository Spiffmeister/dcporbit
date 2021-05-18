module orbit

    using LinearAlgebra
    # using Roots
    using Plots
    using LaTeXStrings
    using JLD2
    using Distributed
    pyplot()

    include("Particle.jl")
    include("OrbitEqns.jl")
    include("Integrators.jl")
    include("Plotting.jl")

    export sim,particle,exact_particle
    export exactsolve,integrate!
    export MagneticForce

end