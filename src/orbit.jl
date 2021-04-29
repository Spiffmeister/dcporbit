module orbit

    using LinearAlgebra
    # using Roots
    using Plots
    using JLD2
    using Distributed
    pyplot()

    include("Particle.jl")
    include("OrbitEqns.jl")
    include("Integrators.jl")

    export exactsolve,integrate!,particle
    export MagneticForce

end