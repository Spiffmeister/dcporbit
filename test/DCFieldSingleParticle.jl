#=
    Test changing fields on a single particle
=#
using dcporbit

mode = FullOrbit

x₀ = [0., 0., 0.01]
v₀ = [1.,0.,0.]

t_f = 2.5
Δt = 1.e-3

α = π/6
B₁(x,t) = [cos(α),-sin(α),0.]
B₂(x,t) = [cos(α),sin(α),0.]

Bfield = [B₁,B₂]

function event(xv1,xv2)
    chk = sign(xv1[3])*sign(xv2[3])
    return chk
end

p = particle(x₀,v₀,FullOrbit,Δt,Bfield,lvol=1,gc_initial=false)
# integrate(p,t_f,eventfn=eventfn)
ODE = forces(Bfield,event=event)
# Run the simulation
solve_orbit!(p,ODE,t_f)
# Plot the orbit
dcporbit.plt_orbit(p)

# mins = findall(x->abs(x)<1.e-15,p.x[3,:])
# @assert length(mins) == 2


pe = analytic_solve(p,Bfield,crossing=true,eventfn=event)


using Plots
pyplot()
plot3d(p.x[1,:],p.x[2,:],p.x[3,:])
plot3d!(pe.x[1,:],pe.x[2,:],pe.x[3,:])



using LinearAlgebra
pn = [norm(cross(p.v[:,i],p.B[:,i]))/norm(p.B[:,i]) for i in 1:length(p.t)]
plot(p.t,pn)


plot(p.t,p.x[3,:])
plot!(pe.t,pe.x[3,:])
scatter!([1.3801464698208694],[0])
