#=
    Performing convergence study
=#

# Is this a guiding center intialisation?
guidingcentre = false
# Set initial conditions
x₀ = [0.,0.,0.]
v₀ = [1.,0.,0.]
# Set simulation parameters
t₀ = 0.
t_f = 1.
Δt = 1.e-1
tol = 1.e-14
# Set the magnetic field
Bfield(x) = [0.,0.,1.]
# Initialise the particle
p = particle(x₀,v₀,guidingcentre,Δt,MagneticForce,Bfield,1)
# Run the simulation
integrate!(p,t_f)
# Compute exact solution
pe = exactsolve(x₀,v₀,p.t,Bfield,crossing=false)


# Plot solutions
p1 = orbit.plt_orbit(p)
p2 = orbit.plt_orbit(pe)