#=
    A simple example of how to use orbit.jl
=#

# Loading the package
using Pkg
Pkg.activate(".")
using orbit



#=
    Simulation parameters
=#
# Use full orbit mode :guidingcentre is the other option
guidingcenter = :fullorbit
# We will initialise at the guiding centre position, we can also set gc_initial=false and set the initial FO position
gc_initial = true
gc₀ = [0., 0., 0.]
# Set the initial velocity
v₀ = [1.,0.,0.]
# The time step size and final time
Δt = 1.e-3
t_f = 30.

#=
    Magnetic field data
=#
# We are going to use a discontinous field so we need an event function which returns a negative when the particle has crossed
function event(xv)
    chk = sign(xv[3,1])*sign(xv[3,2])
    return chk
end

# Define a vector valued function for the magnetic field in each region then store them in a vector of functions
α = π/6
B₁(x) = [cos(α),sin(α),0.]
B₂(x) = [cos(α),-sin(α),0.]
Bfield = [B₁,B₂]
# The initial field the particle will use
lvol = 1

# Construct the struct for the magnetic field data
ODE = forces(Bfield,event=event)

#=
SINGLE PARTICLE case
=#
# Generate the particle data
f₁ = particle(gc₀,v₀,:guidingcenter,Δt,Bfield[1],lvol)
# Solve the particles orbit, note the ! since solve_orbit! is an iterator and appends the data to f₁
solve_orbit!(f₁,ODE,t_f)


#=
MULTI-PARTICLE case
=#

# Set a bunch of initial positions in a vector
x₀ = [[0.,0.,x] for x in -0.9:0.1:-0.1]
# The number of particles we'll be simulating
nparts = length(x₀)

# Create the simulation data, which contains multiple particle objects
f = sim(nparts,x₀,v₀,guidingcenter,Δt,Bfield,1,gc_initial=false)

# Solve all the particle trajectories
run_sim!(f,ODE,t_f)

# Solve the same set of particles but analytically - note that this is _not_ an iterator
fe = analytic_solve(f,Bfield,crossing=true,eventfn=event)



