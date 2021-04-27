module orbit

    using LinearAlgebra
    # using Roots
    using Plots
    using JLD2
    using Distributed
    pyplot()

    include("Integrators.jl")
    include("OrbitEqns.jl")
    include("Particle.jl")

    export exactsolve,integrate

end