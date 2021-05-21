#=
    Performing convergence study
=#

# Is this a guiding center intialisation?
guidingcentre = false
# Set initial conditions
x₀ = [0.,0.,0.]
v₀ = [1.,0.,0.]
# Set simulation parameters
t_f = 1.
Δt = 1.e-3
# Set the magnetic field
Bfield(x) = [0.,0.,1.]
# Initialise the particle
p = particle(x₀,v₀,guidingcentre,Δt,Bfield,1)
ODE = forces(Bfield)
# Run the simulation
solve_orbit!(p,ODE,t_f)

# Compute exact solution
pe1 = analytic_solve(x₀,v₀,p.t,Bfield,crossing=false)

pe2 = analytic_solve(p,Bfield,crossing=false)



# Plot solutions
p1 = orbit.plt_orbit(p)
p2 = orbit.plt_orbit(pe1)