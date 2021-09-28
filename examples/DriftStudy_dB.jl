#=
    Performing convergence study
=#

# Is this a guiding center intialisation?

using Pkg
Pkg.activate(".")
using orbit



function event(xv)
    chk = sign(xv[3,1])*sign(xv[3,2])
    return chk
end

α = π/6
B₁(x) = [cos(α),sin(α),0.]      #0>z
B₂(x) = [cos(α),-sin(α),0.]     #0<z
Bfield = [B₁,B₂]

ODE = forces(Bfield,event=event)


#=
    Guiding center and full orbits moving through Z=0 plane
=#
mode = :fo
Δt = 1.e-3
t_f = 30.

nparts = 4
x₀ = [[0.,0.,-2.],
    [0.,0.,0.],
    [0.,0.,0.],
    [0.,0.,2.]]

v₀ = [1.,0.,0.]

lᵥ = [1,1,2,2]


# Initialise the simulation
f = sim(nparts,x₀,v₀,mode,Δt,Bfield,lᵥ,gc_initial=true)
# Run the simulation
run_sim!(f,ODE,t_f)
# Run the analytic simulation
fe = analytic_solve(f,Bfield,crossing=true,eventfn=event)




#=
Compute the GC positions to add to plot
=#
gc₀ = [0.,0.,0.]
v₀ = [1.,0.,0.]
# Construct GCs
ODE_GC₁ = forces(B₁,force=MagneticForce_GC)
gcsim₁ = particle(gc₀,v₀,:gc,Δt,Bfield[1],1)
ODE_GC₂ = forces(B₂,force=MagneticForce_GC)
gcsim₂ = particle(gc₀,v₀,:gc,Δt,Bfield[2],1)
# Solve the GCs
solve_orbit!(gcsim₁,ODE_GC₁,t_f)
solve_orbit!(gcsim₂,ODE_GC₂,t_f)





#=
=#

using Plots
using Printf
pyplot()


function plt_gcprojection(f,gc₁,gc₂)
    # GC PROJECTION
    cp = palette(:tab20)

    plt = plot(xlabel="x",ylabel="y",legend=true)
    for k = 1:f.nparts
        plot!(plt,f.sp[k].gc[1,:],f.sp[k].gc[2,:],color=cp[k],label="gc_z=$(@sprintf("%.1e",f.sp[k].gc[3,1]))")
    end

    plot!(plt,gc₁.x[1,:],gc₁.x[2,:],linestyle=:dash)
    plot!(plt,gc₂.x[1,:],gc₂.x[2,:],linestyle=:dash)

    return plt
end



function plt_gcavproj(fe,gc₁,gc₂;n=1)
    cp = palette(:tab20)

    plt = plot(xlabel="x",ylabel="y",legend=true,dpi=600,framestyle=:classic)


    gcp = zeros(2,length(f.sp[n].lvol))
    gcp .= NaN
    gcm = zeros(2,length(f.sp[n].lvol))
    gcm .= NaN

    z_p = findall(x->x == 1,f.sp[n].lvol)
    z_m = findall(x->x == 2,f.sp[n].lvol)

    gcp[:,z_p] = f.sp[n].gc[1:2,z_p]
    gcm[:,z_m] = f.sp[n].gc[1:2,z_m]

    plot!(plt,gcp[1,:],gcp[2,:],color=cp[n],label=string("Δt=",f.sp[n].Δt[1],"z<0"))
    plot!(plt,gcm[1,:],gcm[2,:],color=cp[n],label=string("Δt=",f.sp[n].Δt[1],"z>0"))



    t = range(fe.sp[n].t_boundary[2],stop=fe.sp[n].t[end],length=100)
    ave = zeros(3,length(t))
    for i = 1:length(t)
        ave[:,i] = fe.sp[n].avetraj[:,2] * t[i]
    end
    # plot!(plt,ave[1,:],ave[2,:],label=string("gc₀=",f.sp[k].gc[3,1]))
    plot!(plt,ave[1,:],ave[2,:],label="gc_z=$(@sprintf("%.1e",f.sp[n].gc[3,1]))")

    plot!(plt,gc₁.x[1,:],gc₁.x[2,:],linestyle=:dash,label="z<0 field")
    plot!(plt,gc₂.x[1,:],gc₂.x[2,:],linestyle=:dash,label="z>0 field")

    return plt
end






function comp_avtraj(p)

    gcp = zeros(2,size(p.x)[2])
    gcp .= NaN
    gcm = zeros(2,size(p.x)[2])
    gcm .= NaN

    z_p = findall(x->x == 1,p.lvol)
    z_m = findall(x->x == 2,p.lvol)

    z_tmp = findall(x->x != 0,diff(p.lvol))


    gcp[:,setdiff(z_p,z_tmp,z_tmp.+1,z_tmp.-1)] = p.gc[1:2,setdiff(z_p,z_tmp,z_tmp.+1,z_tmp.-1)]
    gcm[:,setdiff(z_m,z_tmp,z_tmp.+1,z_tmp.-1)] = p.gc[1:2,setdiff(z_m,z_tmp,z_tmp.+1,z_tmp.-1)]

    return gcp, gcm
end

# Colour palette
cp = palette(:tab20)

# Plot
pav = plot(xlabel="x",ylabel="y",legend=:topright,framestyle=:classic,dpi=600,foreground_color_grid=:lightgray)

# Particle 1
gc1 = plot!(pav,gcsim₁.x[1,:],gcsim₁.x[2,:],color=cp[1],label="R1 Guiding centre")
fo1 = plot!(pav,f.sp[1].x[1,:],f.sp[1].x[2,:],color=cp[1],label="R1 Full orbit")

# Crossing particle 1
gcp, gcm = comp_avtraj(f.sp[2])
cross1 = plot!(pav,gcp[1,:],gcp[2,:],color=cp[3],label="C1")
plot!(pav,gcm[1,:],gcm[2,:],color=cp[3],label="")
crossav1 = plot!(pav,fe.sp[2].avetraj[1,2]*fe.sp[2].t,fe.sp[2].avetraj[2,2]*fe.sp[2].t,color=cp[3],linestyle=:dash,label="C1 average drift")
crossFO1 = plot!(pav,f.sp[2].x[1,:],f.sp[2].x[2,:],color=cp[3],linestyle=:dashdot,linealpha=0.6,label="C1 full orbit")

# Crossing particle 2
gcp, gcm = comp_avtraj(f.sp[3])
cross2 = plot!(pav,gcp[1,:],gcp[2,:],color=cp[5],label="C2")
plot!(pav,gcm[1,:],gcm[2,:],color=cp[5],label="")
crossav2 = plot!(pav,fe.sp[3].avetraj[1,2]*fe.sp[3].t,fe.sp[3].avetraj[2,2]*fe.sp[3].t,color=cp[5],linestyle=:dash,label="C2 average drift")
crossFO2 = plot!(pav,f.sp[3].x[1,:],f.sp[3].x[2,:],color=cp[5],linestyle=:dashdot,linealpha=0.6,label="C2 full orbit")

# Particle 4
gc2 = plot!(pav,gcsim₂.x[1,:],gcsim₂.x[2,:],color=cp[7],label="R2 Guiding Centre")
fc2 = plot!(pav,f.sp[4].x[1,:],f.sp[4].x[2,:],color=cp[7],label="R2 Full orbit")



savefig(pav,"Figures//gcav.eps")
