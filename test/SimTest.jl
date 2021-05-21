#=
    Test multi-particle
=#

# Number of particles
nparts = 8

# No guiding centre calculations
guidingcenter = false
# Initial GC position
gc₀ = [0., 0., 0.]
# Initial velocity
v₀ = [1.,0.,0.]
# Array of time steps
Δt = [2.0^-i for i in 1:nparts]

t_f = 4.
Bfield(x) = [0.,0.,1.]

# Set the magnetic field for the RHS
dvdt = forces(Bfield)
# Generate the particles
f = sim(nparts,gc₀,v₀,guidingcenter,Δt,Bfield,1)
# Run all particles
run_sim!(f,dvdt,t_f)

g = analytic_solve(f,Bfield,crossing=false)

# Plot an orbit
orbit.plt_orbit(f.sp[1])

orbit.plt_orbit(g.sp[1])