#=
    Performing convergence study
=#

# Is this a guiding center intialisation?


tol = 1.e-14


# Simulation parameters
nparts = 8

guidingcenter = false
gc₀ = [0., 0., 0.]
v₀ = [1.,0.,0.]
Δt = [2.0^-i for i in 1:nparts]
t_f = 40.

α = π/6
B₁(x) = [cos(α),sin(α),0.]
B₂(x) = [cos(α),-sin(α),0.]

Bfield = [B₁,B₂]

function event(xv)
    chk = sign(xv[3,1])*sign(xv[3,2])
    return chk
end


ODE = forces(Bfield,event=event)

f = sim(nparts,gc₀,v₀,guidingcenter,Δt,Bfield,1)

run_sim!(f,ODE,t_f)

g = analytic_solve(f,Bfield,crossing=true,eventfn=event)



using Plots
pyplot()
plot3d(f.sp[4].x[1,:],f.sp[4].x[2,:],f.sp[4].x[3,:])
plot3d!(g.sp[4].x[1,:],g.sp[4].x[2,:],g.sp[4].x[3,:])