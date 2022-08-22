#=
    A simple example of how to use orbit.jl
=#

# Loading the package
# using Pkg
# Pkg.activate(".")
push!(LOAD_PATH,".")
using orbit



#=
    Simulation parameters
=#
# Use full orbit mode :guidingcentre is the other option
guidingcenter = :fullorbit
# We will initialise at the guiding centre position, we can also set gc_initial=false and set the initial FO position
gc_initial = false
# Set the initial position and velocity
x₀ = [0.0,0.0,-0.5]
v₀ = [1.,0.,0.]
# The time step size and final time
Δt = 1.e-3
t_f = 1.

#=
    Magnetic field data
=#
# We are going to use a discontinous field so we need an event function which returns -1 when the particle has crossed
function event(xv)
    chk = sign(xv[3,1])*sign(xv[3,2])
    return chk
end

# Define a vector valued function for the magnetic field in each region then store them in a vector of functions
α = π/6
B₁(x,t) = [cos(α),sin(α),0.]
B₂(x,t) = [cos(α),-sin(α),0.]
Bfield = [B₁,B₂]
# The initial field the particle will use
lvol = 1

# Constructs the magnetic field data
ODE = forces(Bfield,event=event)

#=
SINGLE PARTICLE case
=#
# Generate the particle data
f₁ = particle(x₀,v₀,:fo,Δt,Bfield[1],lvol,gc_initial=false)
# Solve the particles orbit, note the ! since solve_orbit! is an iterator and appends the data to f₁
solve_orbit!(f₁,ODE,t_f)

