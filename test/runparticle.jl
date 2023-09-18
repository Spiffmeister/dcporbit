#=
    Performing convergence study
=#

# using Pkg
# Pkg.activate(".")
using dcporbit

# Is this a guiding center intialisation?
mode = FullOrbit
# Set initial conditions
x₀ = [0.,0.,1.1]
v₀ = [1.,0.,0.]
# Set simulation parameters
t_f = 1.0
Δt = 1.e-3
# Set the magnetic field
α = π/6
Bfield(x,t) = [cos(α),sin(α),0.0]
# Initialise the particle
p = particle(x₀,v₀,mode,Δt,Bfield)
ODE = forces(Bfield)
# Run the simulation
# @time solve_orbit!(p,ODE,t_f,method=:forward_eulers)
@time solve_orbit!(p,ODE,t_f,method=:RK4)

# Compute exact solution
pe = analytic_solve(x₀,v₀,p.t,Bfield,crossing=false)


# Plot solutions
using Plots
pyplot()
plot3d(p.x[1,:],p.x[2,:],p.x[3,:])
plot3d!(pe.x[1,:],pe.x[2,:],pe.x[3,:])


