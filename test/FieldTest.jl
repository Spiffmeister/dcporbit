#=
    Test changing fields
=#

# Is this a guiding center intialisation?


# Simulation parameters
nparts = 8

guidingcenter = false
gc₀ = [0., 0., 0.]
v₀ = [1.,0.,0.]
Δt = 1.e-3

t_f = 5.

α = π/6
B₁(x) = [cos(α),sin(α),0.]
B₂(x) = [cos(α),-sin(α),0.]

Bfield = [B₁,B₂]

function event(xv)
    chk = sign(xv[3,1])*sign(xv[3,2])
    return chk
end


p = particle(gc₀,v₀,guidingcenter,Δt,Bfield,1)
# integrate(p,t_f,eventfn=eventfn)
ODE = forces(Bfield,event=event)
# Run the simulation
solve_orbit!(p,ODE,t_f)
# Plot the orbit
orbit.plt_orbit(p)

mins = findall(x->abs(x)<1.e-15,p.x[3,:])
@assert length(mins) == 2


pe = analytic_solve(p,Bfield,crossing=true,eventfn=event)


# pe1 = analytic_solve(gc₀,v₀,p.t,Bfield,crossing=true,eventfn=event)


println(p.x[3,mins])



using Plots
pyplot()

plot3d(p.x[1,:],p.x[2,:],p.x[3,:])
plot3d!(pe.x[1,:],pe.x[2,:],pe.x[3,:])



ind = findall(x->x!=0,diff(p.lvol))

p.x[3,ind]
pe.x[3,ind]


p.x[:,ind[1]+1]
pe.x[:,ind[1]+1]


