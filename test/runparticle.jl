#=
    Performing convergence study
=#

using Pkg
Pkg.activate(".")
using dcporbit

# Is this a guiding center intialisation?
mode = FullOrbit
# Set initial conditions
x₀ = [0.,0.,0.]
v₀ = [1.,0.,0.]
# Set simulation parameters
t_f = 10.0
Δt = 1.e-3
# Set the magnetic field
Bfield(x,t) = [0.,0.,1.]
# Initialise the particle
p = particle(x₀,v₀,mode,Δt,Bfield)
ODE = forces(Bfield)
# Run the simulation
# @time solve_orbit!(p,ODE,t_f,method=:forward_eulers)
@time solve_orbit!(p,ODE,t_f,method=:RK4)

# Compute exact solution
pe1 = analytic_solve(x₀,v₀,p.t,Bfield,crossing=false)

# pe2 = analytic_solve(p,Bfield,crossing=false)

# Plot solutions
p1 = dcporbit.plt_orbit(p)
# p2 = orbit.plt_orbit(pe1)



